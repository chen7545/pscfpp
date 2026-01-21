/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Reduce.h"
#include "VecOp.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/DeviceMemory.h>
#include <pscf/cuda/cudaErrorCheck.h>
#include <pscf/cuda/complex.h>
#include <util/misc/Log.h>

//#include <thrust/reduce.h>
//#include <thrust/device_vector.h>
//#include <cub/device/device_reduce.cuh>
#include <cub/cub.cuh>

#include <complex>
#include <cmath>

// If defined, use reduction functions from the NVIDIA CUB library
// in preference to hand-coded kernels.
#define USE_NVIDIA_CUB

namespace Pscf {
namespace Reduce {

   // Memory used as workspace by CUB reduction functions
   static DeviceMemory reduceSpace_{};

   // Memory used as workspace for vector transformations
   static DeviceMemory transformSpace_{};

   // CUDA kernels:
   // (defined in anonymous namespace, used within this file only)
   namespace {

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Parallel reduction: single-warp functions
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      /*
      * The reduction algorithm can be simplified during the last 6
      * levels of reduction, because these levels of reduction are
      * performed by a single warp. Within a single warp, each thread
      * executes the same instruction at the same time (SIMD execution).
      * Therefore, we don't need the __syncthreads() command between
      * reduction operations. Further, we do not need to evaluate an
      * if-statement to determine which threads should perform the
      * calculation and which should not, since the entire warp will be
      * dedicated to these operations regardless of whether they perform
      * calculations. Therefore, the if-statement would not free any
      * resources for other tasks, so we omit it for speed.
      *
      * We assume here that a single warp contains 32 threads. All
      * CUDA-compatible GPUs currently meet this criterion, but it is
      * possible that someday there will be GPUs with a different warp
      * size. The methods below may break if the warp size is smaller
      * than 32 threads, because the operations would be performed by
      * multiple warps without __syncthreads() commands to keep them
      * synced. Warps larger than 32 threads would still be compatible
      * with these functions, though the functions are not optimized
      * for this case.
      *
      * These are implemented as separate functions, rather than within
      * the kernels above, because they require the sData array to be
      * defined as volatile (meaning the array values may change at any
      * time, so the compiler must access the actual memory location
      * rather than using cached values).
      */

      #ifndef USE_NVIDIA_CUB
      /*
      * Utility to perform summation reduction within a single warp.
      *
      * \param sData  input array to reduce
      * \param tId  thread ID
      */
      __device__ void _warpSum(volatile cudaReal* sData, int tId)
      {
         sData[tId] += sData[tId + 32];
         sData[tId] += sData[tId + 16];
         sData[tId] += sData[tId + 8];
         sData[tId] += sData[tId + 4];
         sData[tId] += sData[tId + 2];
         sData[tId] += sData[tId + 1];
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Utility to perform maximization reduction within a single warp.
      *
      * \param sData  input array to reduce
      * \param tId  thread ID
      */
      __device__ 
      void _warpMax(volatile cudaReal* sData, int tId)
      {
         if (sData[tId + 32] > sData[tId]) sData[tId] = sData[tId + 32];
         if (sData[tId + 16] > sData[tId]) sData[tId] = sData[tId + 16];
         if (sData[tId + 8] > sData[tId]) sData[tId] = sData[tId + 8];
         if (sData[tId + 4] > sData[tId]) sData[tId] = sData[tId + 4];
         if (sData[tId + 2] > sData[tId]) sData[tId] = sData[tId + 2];
         if (sData[tId + 1] > sData[tId]) sData[tId] = sData[tId + 1];
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Utility to perform minimization reduction within a single warp.
      *
      * \param sData  input array to reduce
      * \param tId  thread ID
      */
      __device__ 
      void _warpMin(volatile cudaReal* sData, int tId)
      {
         if (sData[tId + 32] < sData[tId]) sData[tId] = sData[tId + 32];
         if (sData[tId + 16] < sData[tId]) sData[tId] = sData[tId + 16];
         if (sData[tId + 8] < sData[tId]) sData[tId] = sData[tId + 8];
         if (sData[tId + 4] < sData[tId]) sData[tId] = sData[tId + 4];
         if (sData[tId + 2] < sData[tId]) sData[tId] = sData[tId + 2];
         if (sData[tId + 1] < sData[tId]) sData[tId] = sData[tId + 1];
      }
      #endif

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Parallel reduction: full kernels
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      #ifndef USE_NVIDIA_CUB
      /*
      * Compute sum of array elements (GPU kernel).
      *
      * Assumes each warp is 32 threads.
      * Assumes that each block contains at least 64 threads.
      * Assumes that the block size is a power of 2.
      *
      * \param sum  reduced array containing sum from each thread block
      * \param in  input array
      * \param n  number of input array elements
      */
      __global__ void _sum(cudaReal* sum, const cudaReal* in, int n)
      {
         // number of blocks cut in two to avoid inactive initial threads
         int tId = threadIdx.x;
         int bId = blockIdx.x;
         int bDim = blockDim.x;
         int idx = bId * (bDim*2) + tId;

         // Shared memory holding area
         extern __shared__ cudaReal sData[];

         // Global memory load and first operation
         if (idx < n) {
            sData[tId] = in[idx];
            if (idx + bDim < n) {
               sData[tId] += in[idx+bDim];
            }
         } else {
            // idx > n. Set value to 0.0 to fully populate sData without
            // contributing to the sum
            sData[tId] = 0.0;
         }

         // wait for all threads to finish
         __syncthreads();

         // Make reductions across the block of data, each thread handling
         // one reduction of two data points with strided indices before
         // syncing with each other and then making further reductions.
         for (int stride = bDim / 2; stride > 32; stride /= 2) {
            if (tId < stride) {
               sData[tId] += sData[tId+stride];
            }
            __syncthreads();
         }

         // Unwrap last warp (stride == 32)
         if (tId < 32) {
            _warpSum(sData, tId); // defined at bottom of this file
         }

         // Store the output of the threads in this block
         if (tId == 0) {
            sum[bId] = sData[0];
         }
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Get maximum of array elements (GPU kernel).
      *
      * Assumes each warp is 32 threads.
      * Assumes that each block contains at least 64 threads.
      * Assumes that the block size is a power of 2.
      *
      * \param max  reduced array containing the max from each thread block
      * \param in  input array
      * \param n  number of input array elements
      */
      __global__ void _max(cudaReal* max, const cudaReal* in, int n)
      {
         // number of blocks cut in two to avoid inactive initial threads
         int tId = threadIdx.x;
         int bId = blockIdx.x;
         int bDim = blockDim.x;
         int idx = bId * (bDim*2) + tId;

         // Shared memory holding area
         extern __shared__ cudaReal sData[];

         // Global memory load and first operation
         if (idx < n) {
            sData[tId] = in[idx];
            if (idx + bDim < n) {
               cudaReal in1 = in[idx+bDim];
               sData[tId] = (sData[tId] > in1) ? sData[tId] : in1;
            }
         } else {
            // idx > n. Set value to in[idx-n], an earlier value in the
            // array, to fully populate sData without altering the result
            sData[tId] = in[idx-n];
         }

         // wait for all threads to finish
         __syncthreads();

         // Make reductions across the block of data, each thread handling
         // one reduction of two data points with strided indices before
         // syncing with each other and then making further reductions.
         for (int stride = bDim/2; stride > 32; stride/=2) {
            if (tId < stride) {
               if (sData[tId+stride] > sData[tId]) {
                  sData[tId] = sData[tId+stride];
               }
            }
            __syncthreads();
         }

         // Unwrap last warp (stride == 32)
         if (tId < 32) {
            _warpMax(sData, tId); // defined at bottom of this file
         }

         // Store the output of the threads in this block
         if (tId == 0) {
            max[bId] = sData[0];
         }
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Get maximum absolute magnitude of array elements (GPU kernel).
      *
      * Assumes each warp is 32 threads.
      * Assumes that each block contains at least 64 threads.
      * Assumes that the block size is a power of 2.
      *
      * \param max  reduced array containing the max from each thread block
      * \param in  input array
      * \param n  number of input array elements
      */
      __global__ void _maxAbs(cudaReal* max, const cudaReal* in, int n)
      {
         // number of blocks cut in two to avoid inactive initial threads
         int tId = threadIdx.x;
         int bId = blockIdx.x;
         int bDim = blockDim.x;
         int idx = bId * (bDim*2) + tId;

         // Shared memory holding area
         extern __shared__ cudaReal sData[];

         // Global memory load and first operation
         if (idx < n) {
            sData[tId] = fabs(in[idx]);
            if (idx + bDim < n) {
               cudaReal in1 = fabs(in[idx+bDim]);
               sData[tId] = (sData[tId] > in1) ? sData[tId] : in1;
            }
         } else {
            // idx > n. Set value to 0.0 to fully populate sData without
            // altering the result
            sData[tId] = 0.0;
         }

         // wait for all threads to finish
         __syncthreads();

         // Make reductions across the block of data, each thread handling
         // one reduction across two data points with strided indices before
         // syncing with each other and then making further reductions.
         for (int stride = bDim/2; stride > 32; stride/=2) {
            if (tId < stride) {
               if (sData[tId+stride] > sData[tId]) {
                  sData[tId] = sData[tId+stride];
               }
            }
            __syncthreads();
         }

         // Unwrap last warp (stride == 32)
         if (tId < 32) {
            _warpMax(sData, tId); // defined at bottom of this file
         }

         // Store the output of the threads in this block
         if (tId == 0) {
            max[bId] = sData[0];
         }
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Get minimum of array elements (GPU kernel).
      *
      * Assumes each warp is 32 threads.
      * Assumes that each block contains at least 64 threads.
      * Assumes that the block size is a power of 2.
      *
      * \param min  reduced array containing the min from each thread block
      * \param in  input array
      * \param n  number of input array elements
      */
      __global__ void _min(cudaReal* min, const cudaReal* in, int n)
      {
         // number of blocks cut in two to avoid inactive initial threads
         int tId = threadIdx.x;
         int bId = blockIdx.x;
         int bDim = blockDim.x;
         int idx = bId * (bDim*2) + tId;

         // Shared memory holding area
         extern __shared__ cudaReal sData[];

         // Global memory load and first operation
         if (idx < n) {
            sData[tId] = in[idx];
            if (idx + bDim < n) {
               cudaReal in1 = in[idx+bDim];
               sData[tId] = (sData[tId] < in1) ? sData[tId] : in1;
            }
         } else {
            // idx > n. Set value to in[idx-n], an earlier value in the
            // array, to fully populate sData without altering the result
            sData[tId] = in[idx-n];
         }

         // wait for all threads to finish
         __syncthreads();

         // Make reductions across the block of data, each thread handling
         // one reduction across two data points with strided indices before
         // syncing with each other and then making further reductions.
         for (int stride = bDim/2; stride > 32; stride/=2) {
            if (tId < stride) {
               if (sData[tId+stride] < sData[tId]) {
                  sData[tId] = sData[tId+stride];
               }
            }
            __syncthreads();
         }

         // Unwrap last warp (stride == 32)
         if (tId < 32) {
            _warpMin(sData, tId); // defined at bottom of this file
         }

         // Store the output of the threads in this block
         if (tId == 0) {
            min[bId] = sData[0];
         }
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Get minimum absolute magnitude of array elements (GPU kernel).
      *
      * Assumes each warp is 32 threads.
      * Assumes that each block contains at least 64 threads.
      * Assumes that the block size is a power of 2.
      *
      * \param min  reduced array containing the min from each thread block
      * \param in  input array
      * \param n  number of input array elements
      */
      __global__ void _minAbs(cudaReal* min, const cudaReal* in, int n)
      {
         // number of blocks cut in two to avoid inactive initial threads
         int tId = threadIdx.x;
         int bId = blockIdx.x;
         int bDim = blockDim.x;
         int idx = bId * (bDim*2) + tId;

         // Shared memory holding area
         extern __shared__ cudaReal sData[];

         // Global memory load and first operation
         if (idx < n) {
            sData[tId] = fabs(in[idx]);
            if (idx + bDim < n) {
               cudaReal in1 = fabs(in[idx+bDim]);
               sData[tId] = (sData[tId] < in1) ? sData[tId] : in1;
            }
         } else {
            // idx > n. Set value to fabs(in[idx-n]), an earlier value in the
            // array, to fully populate sData without altering the result
            sData[tId] = fabs(in[idx-n]);
         }

         // wait for all threads to finish
         __syncthreads();

         // Make reductions across the block of data, each thread handling
         // one reduction across two data points with strided indices before
         // syncing with each other and then making further reductions.
         for (int stride = bDim/2; stride > 32; stride/=2) {
            if (tId < stride) {
               if (sData[tId+stride] < sData[tId]) {
                  sData[tId] = sData[tId+stride];
               }
            }
            __syncthreads();
         }

         // Unwrap last warp (stride == 32)
         if (tId < 32) {
            _warpMin(sData, tId); // defined at bottom of this file
         }

         // Store the output of the threads in this block
         if (tId == 0) {
            min[bId] = sData[0];
         }
      }
      #endif

      #ifndef USE_NVIDIA_CUB
      /*
      * Compute inner product of two real arrays (GPU kernel).
      *
      * Assumes each warp is 32 threads.
      * Assumes that each block contains at least 64 threads.
      * Assumes that the block size is a power of 2.
      *
      * \param ip  reduced array with inner prod from each thread block
      * \param a  first input array
      * \param b  second input array
      * \param n  number of input array elements
      */
      __global__ 
      void _innerProduct(cudaReal* ip, 
		         cudaReal const * a,
                         cudaReal const * b, int n)
      {
         // number of blocks cut in two to avoid inactive initial threads
         int tId = threadIdx.x;
         int bId = blockIdx.x;
         int bDim = blockDim.x;
         int idx = bId * (bDim*2) + tId;

         // Shared memory holding area
         extern __shared__ cudaReal sData[];

         // Global memory load and first operation
         if (idx < n) {
            sData[tId] = a[idx] * b[idx];
            if (idx + bDim < n) {
               sData[tId] += (a[idx+bDim] * b[idx+bDim]);
            }
         } else {
            // idx > n. Set value to 0.0 to fully populate sData without
            // contributing to the sum
            sData[tId] = 0.0;
         }

         // wait for all threads to finish
         __syncthreads();

         // Make reductions across the block of data, each thread handling
         // one reduction across two data points with strided indices before
         // syncing with each other and then making further reductions.
         for (int stride = bDim / 2; stride > 32; stride /= 2) {
            if (tId < stride) {
               sData[tId] += sData[tId+stride];
            }
            __syncthreads();
         }

         // Unwrap last warp (stride == 32)
         if (tId < 32) {
            _warpSum(sData, tId); // defined at bottom of this file
         }

         // Store the output of the threads in this block
         if (tId == 0) {
            ip[bId] = sData[0];
         }
      }
      #endif

      // Functors for use in CUB library functions

      #if 0
      // Absolute value functor for a real number.
      struct absRealFunctor {

         __host__ __device__ inline
         cudaReal operator()(cudaReal const & a) const
         {  return std::abs(a); }

      };

      // Binary maximum functor for real numbers.
      struct maxFunctor {

         __host__ __device__ inline
         cudaReal operator() (cudaReal const & a,
                              cudaReal const & b) const
         {  return ((a > b) ? a : b); }

      };

      // Binary minimum functor for real numbers.
      struct minFunctor {

         __host__ __device__ inline
         cudaReal operator() (cudaReal const & a,
                              cudaReal const & b) const
         {  return ((a < b) ? a : b); }

      };

      // Binary addition functor for real numbers.
      struct addRealFunctor {

         __host__ __device__ inline
         cudaReal operator() (cudaReal const & a,
                              cudaReal const & b) const
         {  return a + b; }

      };
      #endif

      // Complex addition functor.
      struct addComplexFunctor {

         __host__ __device__ inline
         cudaComplex operator() (cudaComplex const & a,
                                 cudaComplex const & b) const
         {
            cudaComplex result;
            result.x = a.x + b.x;
            result.y = a.y + b.y;
            return result;
         }

      };

   } // anonmyous namespace Pscf::Reduce::(unnamed)

   // Memory management

   void freeWorkSpace()
   { 
      reduceSpace_.deallocate();
      transformSpace_.deallocate();
   }

   // Public reduction functions

   #ifdef USE_NVIDIA_CUB
   /*
   * Compute sum of elements of a real array.
   *
   * This implementation uses Nvidia CUB Library functions.
   */
   cudaReal sum(DeviceArray<cudaReal> const & a)
   {
      //std::cout << "\n Using CUB implementation of sum";
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      // Create input and output pointers
      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      DeviceArray<cudaReal> out;
      out.allocate(1);
      cudaReal* outPtr = out.cArray();

      // Determine size of required workspace, allocate if needed
      size_t workSize = 0;
      cudaError_t error;
      error = cub::DeviceReduce::Sum(nullptr, workSize,
                                     inPtr, outPtr, n);
      UTIL_CHECK(error == cudaSuccess);
      reduceSpace_.resize(workSize);
      UTIL_CHECK(reduceSpace_.capacity() >= workSize);

      // Perform reduction
      error = cub::DeviceReduce::Sum(reduceSpace_.cArray(), workSize,
                                     inPtr, outPtr, n);
      UTIL_CHECK(error == cudaSuccess);

      // Copy to host and return value
      HostDArray<cudaReal> out_h;
      out_h.allocate(1);
      out_h = out;
      return out_h[0];
   }
   #endif

   #ifndef USE_NVIDIA_CUB
   /*
   * Compute sum of array elements.
   *
   * This implementation uses hand-coded _max kernel.
   */
   cudaReal sum(DeviceArray<cudaReal> const & in)
   {
      //std::cout << "\n Using custom implementation of sum";
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary device arrays for storing reduced data
      DeviceArray<cudaReal> temp1, temp2;

      int i = 0;

      // Perform parallel reduction on GPU repeatedly until n < 1e5
      while (n >= 1e5) {
         // Establish GPU resources for this parallel reduction. Divided by
         // two because of the global memory load in the kernel performing
         // the first level of reduction!
         int nBlocks, nThreads;
         int halvedSize = ceil((float)n/2);

         ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);
         // Note: setThreadsLogical ensures that nThreads is a power of 2

         if (nThreads < 64) {
            // Thread blocks too small. Manually set nThreads to 64
            ThreadArray::setThreadsPerBlock(64);
            ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);

            // If the above was successful, print warning
            Log::file() << "\n Warning: "
                        << "nThreads too small for parallel reduction."
                        << "\n Setting nThreads equal to 64." 
                        << std::endl;
         }

         // Warp size must be 32
         UTIL_CHECK(ThreadArray::warpSize() == 32);

         // Perform parallel reduction
         if (i == 0) { // first reduction, use input array
            temp1.allocate(nBlocks);
            _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), in.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
         } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
            temp2.allocate(nBlocks);
            _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp2.cArray(), temp1.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp1.deallocate();
         } else {                // i is even: reduce temp2, store in temp1
            temp1.allocate(nBlocks);
            _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), temp2.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp2.deallocate();
         }

         n = nBlocks;
         i += 1;
      }

      // Transfer the partially reduced sum to the host
      HostDArray<cudaReal> temp_h;
      if (i == 0) {
         temp_h = in;
      } else if (i % 2 == 1) {
         temp_h = temp1;
      } else {
         temp_h = temp2;
      }

      if (n == 1) {
         return temp_h[0];
      } else {
         // Sum up elements of temp_h to get final result.
         // Use Kahan summation to reduce accumulation of error
         cudaReal sum = 0.0, tempVal, tempSum;
         cudaReal err = 0.0;
         for (int i = 0; i < n; ++i) {
            tempVal = temp_h[i] - err;
            tempSum = sum + tempVal;
            err = tempSum - sum - tempVal;
            sum = tempSum;
         }
         return sum;
      }
   }
   #endif

   /*
   * Compute sum of elements of a complex array.
   *
   * Implementation uses Nvidia CUB Library functions.
   */
   std::complex<cudaReal> sum(DeviceArray<cudaComplex> const & a)
   {
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      // Define pointers to input and output arrays
      cudaComplex* inPtr = const_cast<cudaComplex*>( a.cArray() );
      DeviceArray<cudaComplex> out(1);
      cudaComplex* outPtr = out.cArray();

      // Determine size of required workspace and allocate if necessary
      size_t workSize = 0;
      cudaError_t error;
      auto op = addComplexFunctor{};
      cudaComplex init = makeComplex(0.0, 0.0);
      error = cub::DeviceReduce::Reduce(nullptr, workSize,
                                        inPtr, outPtr, n, op, init);
      UTIL_CHECK(error == cudaSuccess);
      reduceSpace_.resize(workSize);

      // Perform reduction
      error = cub::DeviceReduce::Reduce(reduceSpace_.cArray(), workSize,
                                        inPtr, outPtr, n, op, init);
      UTIL_CHECK(error == cudaSuccess);

      // Copy to host and return value
      HostDArray<cudaComplex> out_h(1);
      out_h = out;
      return std::complex<cudaReal>(out_h[0].x, out_h[0].y);
   }

   /*
   * Return the sum of squares of elements of a real array.
   */
   cudaReal sumSq(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of element-wise squares
      VecOp::sqV(temp, in);

      cudaReal result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }

   /*
   * Return the sum of squares of elements of a complex array.
   */
   std::complex<cudaReal> sumSq(DeviceArray<cudaComplex> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaComplex);
      transformSpace_.resize(workSize);
      DeviceArray<cudaComplex> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of element-wise squares
      VecOp::sqV(temp, in);

      // Compute sum of array temp of element-wise squares
      std::complex<cudaReal> result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }

   #ifdef USE_NVIDIA_CUB
   /*
   * Compute inner product of two real arrays.
   */
   cudaReal innerProduct(DeviceArray<cudaReal> const & a,
                         DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(b.isAllocated());
      UTIL_CHECK(a.capacity() == b.capacity());
      int n = a.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Perform element-wise multiplication v[i] = a[i] * b[i]
      VecOp::mulVV(temp, a, b);

      // Compute sum of element-wise products
      cudaReal result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }
   #endif

   #ifndef USE_NVIDIA_CUB
   /*
   * Compute inner product of two real arrays.
   */
   cudaReal innerProduct(DeviceArray<cudaReal> const & a,
                         DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(b.isAllocated());
      UTIL_CHECK(a.capacity() == b.capacity());
      int n = a.capacity();

      // Set up temporary device arrays for storing reduced data
      DeviceArray<cudaReal> temp1, temp2;

      int i = 0;

      // Perform parallel reduction on GPU repeatedly until n < 1e5
      // (note: for this wrapper, we always call the kernel at least once,
      // even if n < 1e5, so that the part done on the CPU is always just
      // adding up array elements.)
      while (n >= 1e5 || i == 0) {
         // Establish GPU resources for this parallel reduction. Divided by
         // two because of the global memory load in the kernel performing
         // the first level of reduction!
         int nBlocks, nThreads;
         int halvedSize = ceil((float)n/2);

         ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);
         // Note: setThreadsLogical ensures that nThreads is a power of 2

         if (nThreads < 64) {
            // Thread blocks too small. Manually set nThreads to 64
            ThreadArray::setThreadsPerBlock(64);
            ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);

            // If the above was successful, print warning
            Log::file() << "\n Warning: "
                        << "nThreads too small for parallel reduction."
                        << "\n Setting nThreads equal to 64."
			<< std::endl;
         }

         // Warp size must be 32
         UTIL_CHECK(ThreadArray::warpSize() == 32);

         // Perform parallel reduction

         // Note: Only the first kernel call uses kernel _innerProduct.
         // After that, we use the _sum kernel to reduce the array
         // produced by the _innerProduct kernel.

         if (i == 0) { // first reduction, use input arrays
            temp1.allocate(nBlocks);
            _innerProduct<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                               (temp1.cArray(), a.cArray(), b.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
         } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
            temp2.allocate(nBlocks);
            _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp2.cArray(), temp1.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp1.deallocate();
         } else {                // i is even: reduce temp2, store in temp1
            temp1.allocate(nBlocks);
            _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), temp2.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp2.deallocate();
         }

         n = nBlocks;
         i += 1;
      }

      // Transfer the partially reduced sum to the host
      HostDArray<cudaReal> temp_h;
      if (i % 2 == 1) {
         temp_h = temp1;
      } else {
         temp_h = temp2;
      }

      if (n == 1) {
         return temp_h[0];
      } else {
         // Sum up elements of temp_h to get final result.
         // Use Kahan summation to reduce accumulation of error
         cudaReal sum = 0.0, tempVal, tempSum;
         cudaReal err = 0.0;
         for (int i = 0; i < n; ++i) {
            tempVal = temp_h[i] - err;
            tempSum = sum + tempVal;
            err = tempSum - sum - tempVal;
            sum = tempSum;
         }
         return sum;
      }
   }
   #endif

   #ifdef USE_NVIDIA_CUB
   /*
   * Compute max of elements of a real array.
   *
   * This implementation uses an Nvidia CUB library function.
   */
   cudaReal max(DeviceArray<cudaReal> const & a)
   {
      //std::cout << "\n Using CUB implementation of max";
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      DeviceArray<cudaReal> out(1);
      cudaReal* outPtr = out.cArray();

      // Determine size of required workspace, allocated if needed
      size_t workSize = 0;
      cudaError_t error;
      error = cub::DeviceReduce::Max(nullptr, workSize,
                                     inPtr, outPtr, n);
      UTIL_CHECK(error == cudaSuccess);
      reduceSpace_.resize(workSize);

      // Perform reduction
      error = cub::DeviceReduce::Max(reduceSpace_.cArray(), workSize,
                                     inPtr, outPtr, n);
      UTIL_CHECK(error == cudaSuccess);

      // Copy to host and return value
      HostDArray<cudaReal> out_h(1);
      out_h = out;
      return out_h[0];
   }
   #endif

   #ifndef USE_NVIDIA_CUB
   /*
   * Get maximum of array elements.
   *
   * These implementation uses a hand-coded _max kernel.
   */
   cudaReal max(DeviceArray<cudaReal> const & in)
   {
      //std::cout << "\n Using custom implementation of max";
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary device arrays for storing reduced data
      DeviceArray<cudaReal> temp1, temp2;

      int i = 0;

      // Perform parallel reduction on GPU repeatedly until n < 1e5
      while (n >= 1e5) {
         // Establish GPU resources for this parallel reduction. Divided by
         // two because of the global memory load in the kernel performing
         // the first level of reduction!
         int nBlocks, nThreads;
         int halvedSize = ceil((float)n/2);

         ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);
         // Note: setThreadsLogical ensures that nThreads is a power of 2

         if (nThreads < 64) {
            // Thread blocks too small. Manually set nThreads to 64
            ThreadArray::setThreadsPerBlock(64);
            ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);

            // If the above was successful, print warning
            Log::file() << "\n Warning: "
                        << "nThreads too small for parallel reduction."
                        << "\n Setting nThreads equal to 64." 
                        << std::endl;
         }

         // Warp size must be 32
         UTIL_CHECK(ThreadArray::warpSize() == 32);

         // Perform parallel reduction
         if (i == 0) { // first reduction, use input array
            temp1.allocate(nBlocks);
            _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), in.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
         } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
            temp2.allocate(nBlocks);
            _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp2.cArray(), temp1.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp1.deallocate();
         } else {                 // i is even: reduce temp2, store in temp1
            temp1.allocate(nBlocks);
            _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), temp2.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp2.deallocate();
         }

         n = nBlocks;
         i += 1;
      }

      // Transfer the partially reduced sum to the host
      HostDArray<cudaReal> temp_h;
      if (i == 0) {
         temp_h = in;
      } else if (i % 2 == 1) {
         temp_h = temp1;
      } else {
         temp_h = temp2;
      }

      cudaReal max = temp_h[0];
      for (int i = 1; i < n; i++) {
         if (temp_h[i] > max) max = temp_h[i];
      }
      return max;
   }
   #endif

   #ifdef USE_NVIDIA_CUB
   /*
   * Return the maximum absolute magnitude of array elements.
   */
   cudaReal maxAbs(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of absolute magnitudes
      VecOp::absV(temp, in);

      // Compute maximum of array temp of absolute magnitudes
      cudaReal result = Reduce::max(temp);
      temp.dissociate();

      return result;
   }
   #endif

   #ifndef USE_NVIDIA_CUB
   /*
   * Return the maximum absolute magnitude of array elements.
   */
   cudaReal maxAbs(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary device arrays for storing reduced data
      DeviceArray<cudaReal> temp1, temp2;

      int i = 0;

      // Perform parallel reduction on GPU repeatedly until n < 1e5
      //   (note: for this wrapper, we always call the kernel at least once,
      //    even if n < 1e5, so that the part done on the CPU is always just
      //    comparing array element size, without needing fabs().)
      while (n >= 1e5 || i == 0) {
         // Establish GPU resources for this parallel reduction. Divided by
         // two because of the global memory load in the kernel performing
         // the first level of reduction!
         int nBlocks, nThreads;
         int halvedSize = ceil((float)n/2);

         ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);
         // Note: setThreadsLogical ensures that nThreads is a power of 2

         if (nThreads < 64) {
            // Thread blocks too small. Manually set nThreads to 64
            ThreadArray::setThreadsPerBlock(64);
            ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);

            // If the above was successful, print warning
            Log::file() << "\n Warning: "
                        << "nThreads too small for parallel reduction."
                        << "\n Setting nThreads equal to 64." 
                        << std::endl;
         }

         // Warp size must be 32
         UTIL_CHECK(ThreadArray::warpSize() == 32);

         // Perform parallel reduction
         if (i == 0) { // first reduction, use input array
            temp1.allocate(nBlocks);
            _maxAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), in.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
         } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
            temp2.allocate(nBlocks);
            _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp2.cArray(), temp1.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp1.deallocate();
         } else {                 // i is even: reduce temp2, store in temp1
            temp1.allocate(nBlocks);
            _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), temp2.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp2.deallocate();
         }

         n = nBlocks;
         i += 1;
      }

      // Transfer the partially reduced sum to the host
      HostDArray<cudaReal> temp_h;
      if (i == 0) {
         temp_h = in;
      } else if (i % 2 == 1) {
         temp_h = temp1;
      } else {
         temp_h = temp2;
      }

      cudaReal max = temp_h[0];
      for (int i = 1; i < n; i++) {
         if (temp_h[i] > max) max = temp_h[i];
      }
      return max;
   }
   #endif

   #ifdef USE_NVIDIA_CUB
   /*
   * Compute min of elements of a real array.
   *
   * This implementation uses an Nvidia CUB library function.
   */
   cudaReal min(DeviceArray<cudaReal> const & a)
   {
      //std::cout << "\n Using CUB implementation of min";
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      DeviceArray<cudaReal> out;
      out.allocate(1);
      cudaReal* outPtr = out.cArray();

      // Determine size of required workspace, allocated if needed
      size_t workSize = 0;
      cudaError_t error;
      error = cub::DeviceReduce::Min(nullptr, workSize,
                                     inPtr, outPtr, n);
      UTIL_CHECK(error == cudaSuccess);
      reduceSpace_.resize(workSize);

      // Perform reduction
      error = cub::DeviceReduce::Min(reduceSpace_.cArray(), workSize,
                                     inPtr, outPtr, n);
      UTIL_CHECK(error == cudaSuccess);

      // Copy to host and return value
      HostDArray<cudaReal> out_h;
      out_h.allocate(1);
      out_h = out;
      return out_h[0];
   }
   #endif

   #ifndef USE_NVIDIA_CUB
   /*
   * Get minimum of array elements.
   *
   * This implementation uses custom kernels.
   */
   cudaReal min(DeviceArray<cudaReal> const & in)
   {
      //std::cout << "\n Using custom implementation of min";
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary device arrays for storing reduced data
      DeviceArray<cudaReal> temp1, temp2;

      int i = 0;

      // Perform parallel reduction on GPU repeatedly until n < 1e5
      while (n >= 1e5) {
         // Establish GPU resources for this reduction. Divided by two
         // because of the global memory load in the kernel performing
         // the first level of reduction!
         int nBlocks, nThreads;
         int halvedSize = ceil((float)n/2);

         ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);
         // Note: setThreadsLogical ensures that nThreads is a power of 2

         if (nThreads < 64) {
            // Thread blocks too small. Manually set nThreads to 64
            ThreadArray::setThreadsPerBlock(64);
            ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);

            // If the above was successful, print warning
            Log::file() << "\n Warning: "
                        << "nThreads too small for parallel reduction."
                        << "\n Setting nThreads equal to 64." 
                        << std::endl;
         }

         // Warp size must be 32
         UTIL_CHECK(ThreadArray::warpSize() == 32);

         // Perform parallel reduction
         if (i == 0) { // first reduction, use input array
            temp1.allocate(nBlocks);
            _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), in.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
         } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
            temp2.allocate(nBlocks);
            _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp2.cArray(), temp1.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp1.deallocate();
         } else {                 // i is even: reduce temp2, store in temp1
            temp1.allocate(nBlocks);
            _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), temp2.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp2.deallocate();
         }

         n = nBlocks;
         i += 1;
      }

      // Transfer the partially reduced sum to the host
      HostDArray<cudaReal> temp_h;
      if (i == 0) {
         temp_h = in;
      } else if (i % 2 == 1) {
         temp_h = temp1;
      } else {
         temp_h = temp2;
      }

      cudaReal min = temp_h[0];
      for (int i = 1; i < n; i++) {
         if (temp_h[i] < min) min = temp_h[i];
      }
      return min;
   }
   #endif

   #ifdef USE_NVIDIA_CUB
   /*
   * Return the minimum absolute magnitude of array elements.
   */
   cudaReal minAbs(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of absolute magnitudes
      VecOp::absV(temp, in);

      // Compute minimum of array temp of absolute magnitudes
      cudaReal result = Reduce::min(temp);
      temp.dissociate();

      return result;
   }
   #endif

   #ifndef USE_NVIDIA_CUB
   /*
   * Get minimum absolute magnitude of array elements.
   */
   cudaReal minAbs(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary device arrays for storing reduced data
      DeviceArray<cudaReal> temp1, temp2;

      int i = 0;

      // Perform parallel reduction on GPU repeatedly until n < 1e5
      //   (note: for this wrapper, we always call the kernel at least once,
      //    even if n < 1e5, so that the part done on the CPU is always just
      //    comparing array element size, without needing fabs().)
      while (n >= 1e5 || i == 0) {
         // Establish GPU resources for this parallel reduction. Divided by
         // two because of the global memory load in the kernel performing
         // the first level of reduction!
         int nBlocks, nThreads;
         int halvedSize = ceil((float)n/2);

         ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);
         // Note: setThreadsLogical ensures that nThreads is a power of 2

         if (nThreads < 64) {
            // Thread blocks too small. Manually set nThreads to 64
            ThreadArray::setThreadsPerBlock(64);
            ThreadArray::setThreadsLogical(halvedSize,nBlocks,nThreads);

            // If the above was successful, print warning
            Log::file() << "\n Warning: "
                        << "nThreads too small for parallel reduction."
                        << "\n Setting nThreads equal to 64." 
                        << std::endl;
         }

         // Warp size must be 32
         UTIL_CHECK(ThreadArray::warpSize() == 32);

         // Perform parallel reduction
         if (i == 0) { // first reduction, use input array
            temp1.allocate(nBlocks);
            _minAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), in.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
         } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
            temp2.allocate(nBlocks);
            _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp2.cArray(), temp1.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp1.deallocate();
         } else {                 // i is even: reduce temp2, store in temp1
            temp1.allocate(nBlocks);
            _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                       (temp1.cArray(), temp2.cArray(), n);
            cudaErrorCheck( cudaGetLastError() );
            temp2.deallocate();
         }

         n = nBlocks;
         i += 1;
      }

      // Transfer the partially reduced sum to the host
      HostDArray<cudaReal> temp_h;
      if (i == 0) {
         temp_h = in;
      } else if (i % 2 == 1) {
         temp_h = temp1;
      } else {
         temp_h = temp2;
      }

      cudaReal min = temp_h[0];
      for (int i = 1; i < n; i++) {
         if (temp_h[i] < min) min = temp_h[i];
      }
      return min;
   }
   #endif

} // namespace Reduce
} // namespace Pscf

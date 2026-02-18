/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOpFts.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rpg {
namespace VecOpFts {

   // CUDA kernels: 
   // (defined in anonymous namespace, used only in this file)

   namespace {

      // Rescale array a from [0,1] to [-b, b]
      __global__ void _mcftsScale(cudaReal* a, cudaReal const b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = a[i] * 2 * b - b;
         }
      }

      // Add array b to real part of a and array c to imaginary part of a
      __global__ void _fourierMove(cudaComplex* a, cudaReal const * b, 
                                   cudaReal const * c, const int n) {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i];
            a[i].y += c[i];
         }
      }

      // Compute d field (functional derivative of H[w])
      __global__ void _computeDField(cudaReal* d, cudaReal const * Wc, 
                                     cudaReal const * Cc, cudaReal const a, 
                                     cudaReal const b, cudaReal const s, 
                                     const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            d[i] = a * (b * (Wc[i] - s) + Cc[i]);
         }
      }

      // Compute force bias
      __global__ void _computeForceBias(cudaReal* result, cudaReal const * di, 
                                        cudaReal const * df, 
                                        cudaReal const * dwc, 
                                        cudaReal mobility, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            result[i] = 0.5 * (di[i] + df[i]) * 
                        (dwc[i] + mobility * (0.5 * (di[i] - df[i])));
         }
      }

      // Shift w Field
      template <int D>
      __global__ void _shiftWField(cudaReal*  wshift,
                                   cudaReal const * w0,
                                   int const * meshDims, 
                                   int const * shift,
                                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         int position[D];
         int shiftPosition[D];
         int offsets[D];
         offsets[D - 1] = 1;
         int shiftRank;
         for (int d = D - 1; d > 0; --d) {
            offsets[d-1] = offsets[d]*meshDims[d];
         }

         for (int i = startID; i < n; i += nThreads) {

            // Convert index i to position
            int remainder = i;
            int j;
            for (j = 0; j < D - 1; ++j){
               position[j] = remainder/offsets[j];
               remainder -= position[j]*offsets[j];
            }
            position[j] = remainder;

            for (int d = 0; d < D; ++d){
               if (shift[d]>= 0){
                  shiftPosition[d] = (position[d] + shift[d]) % meshDims[d];
               }  else{
                  shiftPosition[d] = (position[d] + shift[d] + meshDims[d]) % meshDims[d];
               }
            }

            // Calculate the rank after shifting
            shiftRank = 0;
            for (j = 0; j < D -1; j++){
               shiftRank += shiftPosition[j]*offsets[j];
            }
            shiftRank += shiftPosition[j];

            wshift[shiftRank] = w0[i];
         }
      }


   }

   // Kernel wrappers:

   /*
   * Rescale array a from [0,1] to [-b, b], GPU kernel wrapper.
   */
   void mcftsScale(DeviceArray<cudaReal>& a, cudaReal const b)
   {
      const int n = a.capacity();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mcftsScale<<<nBlocks, nThreads>>>(a.cArray(), b, n);
   }

   /*
   * Add array b to real part of a and array c to imaginary part of a
   */
   void fourierMove(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _fourierMove<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), 
                                          c.cArray(), n);
   }

   /*
   * Compute d field (functional derivative of H[w])
   */
   void computeDField(DeviceArray<cudaReal>& d, 
                      DeviceArray<cudaReal> const & Wc, 
                      DeviceArray<cudaReal> const & Cc, 
                      cudaReal const a, cudaReal const b, cudaReal const s)
   {
      const int n = d.capacity();
      UTIL_CHECK(Wc.capacity() == n);
      UTIL_CHECK(Cc.capacity() == n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _computeDField<<<nBlocks, nThreads>>>(d.cArray(), Wc.cArray(), 
                                            Cc.cArray(), a, b, s, n);
   }

   /*
   * Compute force bias
   */
   void computeForceBias(DeviceArray<cudaReal>& result, 
                         DeviceArray<cudaReal> const & di, 
                         DeviceArray<cudaReal> const & df, 
                         DeviceArray<cudaReal> const & dwc, cudaReal mobility)
   {
      const int n = result.capacity();
      UTIL_CHECK(di.capacity() == n);
      UTIL_CHECK(df.capacity() == n);
      UTIL_CHECK(dwc.capacity() == n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _computeForceBias<<<nBlocks, nThreads>>>(result.cArray(), di.cArray(), 
                                               df.cArray(), dwc.cArray(), 
                                               mobility, n);
   }

   /**
   * Shift w Field
   */
   template <int D>
   void shiftWField(DeviceArray<cudaReal>& wshift,
                    DeviceArray<cudaReal>const & w0,
                    IntVec<D> const & meshDims,
                    IntVec<D> shift)
   {
      const int n = wshift.capacity();
      UTIL_CHECK(w0.capacity() == n);

      HostDArray<int> meshDims_h(D);
      HostDArray<int> shift_h(D);

      for (int i = 0; i < D; ++i) {
         meshDims_h[i] = meshDims[i];
         shift_h[i] = shift[i];
      }

      DeviceArray<int> meshDims_d(D);
      DeviceArray<int> shift_d(D);
      meshDims_d = meshDims_h;
      shift_d = shift_h;
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _shiftWField<D><<<nBlocks, nThreads>>>(wshift.cArray(), w0.cArray(), 
                                             meshDims_d.cArray(), shift_d.cArray(), n);

   }

   // Instantiate shiftWField methods for D = 1, 2, 3
   template void shiftWField<1>(DeviceArray<cudaReal>& wshift,
                                DeviceArray<cudaReal>const & w0,
                                IntVec<1> const & meshDims,
                                IntVec<1> shift);
   template void shiftWField<2>(DeviceArray<cudaReal>& wshift,
                                DeviceArray<cudaReal>const & w0,
                                IntVec<2> const & meshDims,
                                IntVec<2> shift);
   template void shiftWField<3>(DeviceArray<cudaReal>& wshift,
                                DeviceArray<cudaReal>const & w0,
                                IntVec<3> const & meshDims,
                                IntVec<3> shift);

}
}
}

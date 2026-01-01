#ifndef PSCF_CUDA_REDUCE_H
#define PSCF_CUDA_REDUCE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>

namespace Pscf {

   /**
   * Functions that perform parallel reductions on the GPU.
   *
   * A reduction is any operation that involves reducing all of the
   * elements of an array (or a set of multiple arrays) down to a single
   * result. Examples include taking the sum or finding the maximum of
   * all array elements. Highly efficient algorithms have been developed
   * to perform such operations in parallel on a GPU, and those
   * algorithms are implemented here.
   *
   * A kernel wrapper is provided for each reduction operation, which
   * takes a DeviceArray as an input and returns a single value. The
   * wrappers are called on the host, and they call CUDA kernels
   * internally to perform the reductions in parallel.
   *
   * The CUDA kernels that perform the parallel reductions are not
   * intended to be called directly. The wrappers and CUDA kernels are
   * implemented together in a way that optimizes the speed of the
   * overall operation. To do this, several important aspects of the
   * algorithm are performed by the wrapper, detailed further below.
   *
   * First, the CUDA kernels do not perform the entire reduction.
   * Rather, each block of threads in the GPU is reduced to a single
   * value, and the kernel returns an array containing one value for
   * each thread block. The kernel wrapper then performs the remaining
   * reduction by either calling the kernel again or performing the
   * remaining reduction on the host, depending on the size of the
   * array output by the first kernel. Once the array is smaller than
   * 1e5 elements, the remaining reduction is performed on the CPU.
   * (This value is chosen because the GPU does not significantly 
   * speed up reductions on arrays smaller than 1e5 elements.)
   *
   * Second, the CUDA kernels are designed to use a number of threads
   * that is only half the size of the input array (rounded up to the
   * nearest multiple of the thread block size). The kernel wrapper
   * properly determines the GPU configuration before calling the kernel.
   *
   * \ingroup Pscf_Cuda_Module
   * @{
   */
   namespace Reduce {
  
      // Summation
      
      /**
      * Compute sum of real array elements.
      *
      * \param in  input array
      */
      cudaReal sum(DeviceArray<cudaReal> const & in);
   
      /**
      * Compute sum of array of complex elements.
      *
      * \param in  input array
      */
      cudaComplex sum(DeviceArray<cudaComplex> const & in);
   
      // Vector inner products 
     
      /**
      * Compute inner product of two real arrays.
      *
      * \param a  first input array
      * \param b  second input array
      */
      cudaReal innerProduct(DeviceArray<cudaReal> const & a,
                            DeviceArray<cudaReal> const & b);
   
      // Extrema - real elements
     
      /**
      * Get maximum of array elements.
      *
      * \param in  input array
      */
      cudaReal max(DeviceArray<cudaReal> const & in);
   
      /**
      * Get maximum absolute magnitude of real array elements.
      *
      * \param in  input array
      */
      cudaReal maxAbs(DeviceArray<cudaReal> const & in);
   
      /**
      * Get minimum of array elements.
      *
      * \param in  input array
      */
      cudaReal min(DeviceArray<cudaReal> const & in);
   
      /**
      * Get minimum absolute magnitude of real array elements.
      *
      * \param in  input array
      */
      cudaReal minAbs(DeviceArray<cudaReal> const & in);
   
      // Memory management
   
      /**
      * Free any work space currently allocated for reductions.
      */
      void freeWorkSpace(); 
   
      /** @} */
   
   } // namespace Reduce
} // namespace Pscf
#endif

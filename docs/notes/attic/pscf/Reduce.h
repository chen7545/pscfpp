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
#include <complex>

namespace Pscf {

   /**
   * Functions that perform parallel reductions on a GPU.
   *
   * A reduction is any operation that involves reducing all of the
   * elements of an array (or several arrays) to a single scalar result. 
   * Examples include taking the sum or finding the maximum of all
   * array elements. Highly efficient algorithms have been developed
   * to perform such operations in parallel on a GPU, and those
   * algorithms are implemented here.
   *
   * A wrapper function is provided for each reduction operation, which
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
   * \defgroup Pscf_Cuda_Reduce_Module Reduce (GPU)
   * \ingroup Pscf_Cuda_Module
   */

   /**
   * Reduction operations on CPU or GPU.
   *
   * \ingroup Pscf_Math_Module
   */
   namespace Reduce {

      // Summation

      /**
      * Return sum of elements of a real array.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of elements
      */
      cudaReal sum(DeviceArray<cudaReal> const & in);

      /**
      * Return sum of elements of a complex array.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  complex input array
      * \return  complex sum of elements
      */
      std::complex<cudaReal> sum(DeviceArray<cudaComplex> const & in);

      // Sum of squares and products 

      /**
      * Return sum of squares of elements of a real array.
      *
      * This function returns the square of the Euclidean norm.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of elements
      */
      cudaReal sumSq(DeviceArray<cudaReal> const & in);

      /**
      * Return sum of squares of elements of a complex array.
      *
      * This function returns the complex sum of complex squares of elements.
      * This is not the same as the Hilbert space inner product.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of elements
      */
      std::complex<cudaReal> sumSq(DeviceArray<cudaComplex> const & in);

      /**
      * Return inner product of two real arrays.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param a  first real input array
      * \param b  second real input array
      * \return  Euclidean inner product
      */
      cudaReal innerProduct(DeviceArray<cudaReal> const & a,
                            DeviceArray<cudaReal> const & b);

      // Extrema - real array inputs

      /**
      * Return maximum of real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  value of maximum element
      */
      cudaReal max(DeviceArray<cudaReal> const & in);

      /**
      * Get maximum absolute magnitude of real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  magnitude of element of maximum magnitude
      */
      cudaReal maxAbs(DeviceArray<cudaReal> const & in);

      /**
      * Return minimum of real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  value of minimum element
      */
      cudaReal min(DeviceArray<cudaReal> const & in);

      /**
      * Return minimum absolute magnitude of real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  magnitude of element of minimum magnitude
      */
      cudaReal minAbs(DeviceArray<cudaReal> const & in);

      // Memory management

      /**
      * Free any private work space currently allocated for reductions.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      */
      void freeWorkSpace();

   } // namespace Reduce
} // namespace Pscf
#endif

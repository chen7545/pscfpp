#ifndef PSCF_REDUCE_CX_H
#define PSCF_REDUCE_CX_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>
#include <complex.h>

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Functions that perform reductions on fftw_complex arrays.
   *
   * A reduction is any operation that involves reducing all of 
   * the elements of one or more array to a single scalar, such as
   * a summation or inner product. 
   *
   * \ingroup Prdc_Cpu_Module
   * @{
   */
   namespace Reduce {

      /**
      * Compute sum of array elements (complex).
      *
      * \param a  input array
      * \return complex sum of elements of array a
      */
      std::complex<double> sum(Array<fftw_complex> const & a);

      /**
      * Compute sum of squares of elements of a complex array.
      *
      * \param a  input array
      * \return complex sum of squares of elements of a 
      */
      std::complex<double> sumSq(Array<fftw_complex> const & a);

      /**
      * Compute sum of complex products of elements of two arrays (complex).
      *
      * \param a  first input array
      * \param b  second input array
      * \return sum of elementwise complex products a[i]*b[i]
      */
      std::complex<double> sumProduct(Array<fftw_complex> const & a,
                                      Array<fftw_complex> const & b);

   /** @} */

   } // namespace Reduce
} // namespace Pscf
#endif

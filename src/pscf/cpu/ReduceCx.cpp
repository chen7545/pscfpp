/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ReduceCx.h"
#include <util/containers/Array.h>
#include <fftw3.h>
#include <complex.h>

namespace Pscf {
namespace Reduce {

   /*
   * Compute sum of complex array elements.
   */
   std::complex<double> sum(Array<fftw_complex> const & a)
   {
      int n = a.capacity();
      UTIL_CHECK(n > 0);
      fftw_complex sum;
      sum[0] = 0.0;
      sum[1] = 0.0;
      for (int i = 0; i < n; i++) {
         sum[0] += a[i][0];
         sum[1] += a[i][1];
      }
      std::complex<double> out(sum[0], sum[1]);
      return out;
   }

   /*
   * Compute sum of square of elements of a complex array.
   */
   std::complex<double> sumSq(Array<fftw_complex> const & a)
   {
      int n = a.capacity();
      UTIL_CHECK(n > 0);
      fftw_complex sum;
      double ar, ai;
      sum[0] = 0.0;
      sum[1] = 0.0;
      for (int i = 0; i < n; i++) {
         ar = a[i][0];
         ai = a[i][1];
         sum[0] += ar*ar - ai*ai;
         sum[1] += 2.0 * ar * ai;
      }
      std::complex<double> out(sum[0], sum[1]);
      return out;
   }

   /*
   * Compute sum of element-wise product of two complex arrays.
   */
   std::complex<double> sumProduct(Array<fftw_complex> const & a,
                                   Array<fftw_complex> const & b)
   {
      int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      fftw_complex sum;
      sum[0] = 0.0;
      sum[1] = 0.0;
      for (int i = 0; i < n; i++) {
         sum[0] += a[i][0] * b[i][0] - a[i][1] * b[i][1];
         sum[1] += a[i][0] * b[i][1] + a[i][1] * b[i][0];
      }
      std::complex<double> out(sum[0], sum[1]);
      return out;
   }

} // Reduce
} // Pscf

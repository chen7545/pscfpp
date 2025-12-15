/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOpCx.h"
#include <util/containers/Array.h>
#include <fftw3.h>
#include <cmath>

namespace Pscf {
namespace Prdc {
namespace Cpu {
namespace VecOp {

   using namespace Util;

   // Assignment

   /*
   * Vector-vector assignment, a[i] = b[i] (complex).
   */
   void eqV(Array<fftw_complex> & a, Array<fftw_complex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0];
         a[i][1] = b[i][1];
      }
   }

   /*
   * Vector-vector assignment, a[i] = b[i] (mixed).
   */
   void eqV(Array<fftw_complex> & a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i];
         a[i][1] = 0.0;
      }
   }

   /*
   * Vector-scalar assignment, a[i][0] = b (complex).
   */
   void eqS(Array<fftw_complex> & a, fftw_complex b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[0];
         a[i][1] = b[1];
      }
   }

   /*
   * Vector-scalar assignment, a[i] = b (mixed).
   */
   void eqS(Array<fftw_complex> & a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b;
         a[i][1] = 0.0;
      }
   }

   // Addition

   /*
   * Vector-vector addition, a[i] = b[i] + c[i] (complex).
   */
   void addVV(Array<fftw_complex> & a,
              Array<fftw_complex> const & b, 
              Array<fftw_complex> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] + c[i][0];
         a[i][1] = b[i][1] + c[i][1];
      }
   }

   /*
   * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
   */
   void addVV(Array<fftw_complex> & a,
              Array<fftw_complex> const & b, 
              Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] + c[i];
         a[i][1] = b[i][1];
      }
   }

   /*
   * Vector-scalar addition, a[i][0] = b[i][0] + c (complex).
   */
   void addVS(Array<fftw_complex> & a, 
              Array<fftw_complex> const & b, 
              fftw_complex c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);

      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] + c[0];
         a[i][1] = b[i][1] + c[1];
      }
   }

   /*
   * Vector-scalar addition, a[i][0] = b[i][0] + c (mixed).
   */
   void addVS(Array<fftw_complex> & a, 
              Array<fftw_complex> const & b, 
              double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] + c;
         a[i][1] = b[i][1];
      }
   }

   // Subtraction

   /*
   * Vector-vector subtraction, a[i] = b[i] - c[i] (complex).
   */
   void subVV(Array<fftw_complex> & a,
              Array<fftw_complex> const & b, 
              Array<fftw_complex> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] - c[i][0];
         a[i][1] = b[i][1] - c[i][1];
      }
   }

   /*
   * Vector-vector subtraction, a[i] = b[i] - c[i] (mixed).
   */
   void subVV(Array<fftw_complex> & a,
              Array<fftw_complex> const & b, 
              Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] - c[i];
         a[i][1] = b[i][1];
      }
   }

   /*
   * Vector-scalar subtraction, a[i] = b[i] - c (complex).
   */
   void subVS(Array<fftw_complex>& a,
              Array<fftw_complex> const & b,
              fftw_complex c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] - c[0];
         a[i][1] = b[i][1] - c[1];
      }
   }

   /*
   * Vector-scalar subtraction, a[i] = b[i] - c (mixed).
   */
   void subVS(Array<fftw_complex> & a,
              Array<fftw_complex> const & b,
              double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] - c;
         a[i][1] = b[i][1];
      }
   }

   // Multiplication

   /*
   * Vector-vector multiplication, a[i] = b[i] * c[i] (complex).
   */
   void mulVV(Array<fftw_complex>& a,
              Array<fftw_complex> const & b, 
              Array<fftw_complex> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0]*c[i][0] - b[i][1]*c[i][1];
         a[i][1] = b[i][0]*c[i][1] + b[i][1]*c[i][0];
      }
   }

   /*
   * Vector-vector multiplication, a[i] = b[i] * c[i] (mixed).
   */
   void mulVV(Array<fftw_complex> & a,
              Array<fftw_complex> const & b, 
              Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0]*c[i];
         a[i][1] = b[i][1]*c[i];
      }
   }

   /*
   * Vector-scalar multiplication, a[i] = b[i] * c (complex)
   */
   void mulVS(Array<fftw_complex> & a,
              Array<fftw_complex> const & b,
              fftw_complex c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0]*c[0] - b[i][1]*c[1];
         a[i][1] = b[i][0]*c[1] + b[i][1]*c[0];
      }
   }

   /*
   * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
   */
   void mulVS(Array<fftw_complex>& a,
              Array<fftw_complex> const & b,
              double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0]*c;
         a[i][1] = b[i][1]*c;
      }
   }

   // Division

   /*
   * Vector-vector division, a[i] = b[i] / c[i] (complex).
   */
   void divVV(Array<fftw_complex>& a,
              Array<fftw_complex> const & b, 
              Array<fftw_complex> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      double cr, ci, cSq;
      for (int i = 0; i < n; ++i) {
         cr = c[i][0]; 
         ci = c[i][1]; 
         cSq = cr*cr + ci*ci;
         a[i][0] = (b[i][0]*cr + b[i][1]*ci)/cSq;
         a[i][1] = (b[i][1]*cr - b[i][0]*ci)/cSq;
      }
   }

   /*
   * Vector-vector division, a[i] = b[i] / c[i] (mixed).
   */
   void divVV(Array<fftw_complex>& a,
              Array<fftw_complex> const & b, 
              Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] / c[i];
         a[i][1] = b[i][1] / c[i];
      }
   }

   /*
   * Vector-scalar division, a[i][0] = b[i][0] / c (complex).
   */
   void divVS(Array<fftw_complex>& a, 
                    Array<fftw_complex> const & b, 
                    fftw_complex c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      double cr = c[0];
      double ci = c[1];
      double cSq = cr * cr + ci * ci;
      for (int i = 0; i < n; ++i) {
         a[i][0] = (b[i][0]*cr + b[i][1]*ci)/cSq;
         a[i][1] = (b[i][1]*cr - b[i][0]*ci)/cSq;
      }
   }

   /*
   * Vector-scalar division, a[i][0] = b[i][0] / c (mixed).
   */
   void divVS(Array<fftw_complex>& a, 
                    Array<fftw_complex> const & b, 
                    double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] = b[i][0] / c;
         a[i][1] = b[i][1] / c;
      }
   }

   // Exponentiation

   /*
   * Vector exponentiation, a[i] = exp(b[i]) (complex).
   */
   void expV(Array<fftw_complex>& a, Array<fftw_complex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      double rex;
      for (int i = 0; i < n; ++i) {
         rex = std::exp(b[i][0]);
         a[i][0] = rex * cos(b[i][1]);
         a[i][1] = rex * sin(b[i][1]);
      }
   }

   // Compound addition

   /*
   * Vector-vector in-place addition, a[i][0] += b[i][0] (complex).
   */
   void addEqV(Array<fftw_complex> & a, Array<fftw_complex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] += b[i][0];
         a[i][1] += b[i][1];
      }
   }

   /*
   * Vector-vector in-place addition, a[i][0] += b[i][0] (mixed).
   */
   void addEqV(Array<fftw_complex> & a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] += b[i];
      }
   }

   /*
   * Vector-scalar in-place addition, a[i][0] += b (complex).
   */
   void addEqS(Array<fftw_complex>& a, fftw_complex b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] += b[0];
         a[i][1] += b[1];
      }
   }

   /*
   * Vector-scalar in-place addition, a[i][0] += b (mixed).
   */
   void addEqS(Array<fftw_complex>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] += b;
      }
   }

   // Compound subtraction

   /*
   * Vector-vector in-place subtraction, a[i] -= b[i] (complex).
   */
   void subEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] -= b[i][0];
         a[i][1] -= b[i][1];
      }
   }

   /*
   * Vector-vector in-place subtraction, a[i] -= b[i] (mixed).
   */
   void subEqV(Array<fftw_complex>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] -= b[i];
      }
   }

   /*
   * Vector-scalar in-place subtraction, a[i][0] -= b (complex).
   */
   void subEqS(Array<fftw_complex>& a, fftw_complex b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] -= b[0];
         a[i][1] -= b[1];
      }
   }

   /*
   * Vector-scalar in-place subtraction, a[i][0] -= b (mixed).
   */
   void subEqS(Array<fftw_complex>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] -= b;
      }
   }

   // Compound multiplication

   /*
   * Vector-vector in-place multiplication, a[i] *= b[i] (complex).
   */
   void mulEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      double x, y;
      for (int i = 0; i < n; ++i) {
         x = a[i][0] * b[i][0] - a[i][1] * b[i][1];
         y = a[i][0] * b[i][1] + a[i][1] * b[i][0];
         a[i][0] = x;
         a[i][1] = y;
      }
   }

   /*
   * Vector-vector in-place multiplication, a[i] *= b[i] (mixed).
   */
   void mulEqV(Array<fftw_complex>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] *= b[i];
         a[i][1] *= b[i];
      }
   }

   /*
   * Vector-scalar in-place multiplication, a[i] *= b[i] (complex).
   */
   void mulEqV(Array<fftw_complex>& a, fftw_complex const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      double x, y;
      for (int i = 0; i < n; ++i) {
         x = a[i][0] * b[0] - a[i][1] * b[1];
         y = a[i][0] * b[1] + a[i][1] * b[0];
         a[i][0] = x;
         a[i][1] = y;
      }
   }

   /*
   * Vector-scalar in-place multiplication, a[i][0] *= b (mixed)
   */
   void mulEqS(Array<fftw_complex>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] *= b;
         a[i][1] *= b;
      }
   }

   // Compound division

   /*
   * Vector-vector in-place division, a[i] /= b[i] (complex).
   */
   void divEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      double br, bi, bSq, x, y;
      for (int i = 0; i < n; ++i) {
         br = b[i][0]; 
         bi = b[i][1]; 
         bSq = br*br + bi*bi;
         x = (a[i][0] * br + a[i][1] * bi)/bSq;
         y = (a[i][1] * br - a[i][0] * bi)/bSq;
         a[i][0] = x;
         a[i][1] = y;
      }
   }

   /*
   * Vector-vector in-place division, a[i] /= b[i] (mixed).
   */
   void divEqV(Array<fftw_complex>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i][0] /= b[i];
         a[i][1] /= b[i];
      }
   }

   /*
   * Vector-scalar in-place division, a[i] /= b (complex).
   */
   void divEqS(Array<fftw_complex>& a, fftw_complex b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      double br = b[0]; 
      double bi = b[1]; 
      double bSq = br*br + bi*bi;
      double x, y;
      for (int i = 0; i < n; ++i) {
         x = (a[i][0] * br + a[i][1] * bi) / bSq;
         y = (a[i][1] * br - a[i][0] * bi) / bSq;
         a[i][0] = x;
         a[i][1] = y;
      }
   }

   /*
   * Vector-scalar in-place division, a[i] /= b (mixed).
   */
   void divEqS(Array<fftw_complex>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i][0] /= b;
         a[i][1] /= b;
      }
   }

} // namespace VecOp
} // namespace Cpu
} // namespace Prdc
} // namespace Pscf

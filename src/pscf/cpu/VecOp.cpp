/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "VecOp.h"
#include <util/containers/Array.h>
#include <cmath>

namespace Pscf {
namespace VecOp {

   // Assignment

   /*
   * Vector-vector assignment, a[i] = b[i] (real).
   */
   void eqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i];
      }
   }

   /*
   * Vector-scalar assignment, a[i] = b (real).
   */
   void eqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] = b;
      }
   }

   // Addition

   /*
   * Vector-vector addition, a[i] = b[i] + c[i] (real).
   */
   void addVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] + c[i];
      }
   }

   /*
   * Vector-scalar addition, a[i] = b[i] + c (real).
   */
   void addVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);

      for (int i = 0; i < n; ++i) {
         a[i] = b[i] + c;
      }
   }

   // Subtraction

   /*
   * Vector-vector subtraction, a[i] = b[i] - c[i] (real).
   */
   void subVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] - c[i];
      }
   }

   /*
   * Vector-scalar subtraction, a[i] = b[i] - c (real).
   */
   void subVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] - c;
      }
   }

   // Multiplication

   /*
   * Vector-vector multiplication, a[i] = b[i] * c[i] (real).
   */
   void mulVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] * c[i];
      }
   }

   /*
   * Vector-scalar multiplication, a[i] = b[i] * c (real).
   */
   void mulVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] * c;
      }
   }

   // Division

   /*
   * Vector-vector division, a[i] = b[i] / c[i] (real).
   */
   void divVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] / c[i];
      }
   }

   /*
   * Vector-scalar division, a[i] = b[i] / c (real).
   */
   void divVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] / c;
      }
   }

   /*
   * Scalar-vector division, a[i] = b / c[i] (real).
   */
   void divSV(Array<double>& a, double b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b/c[i];
      }
   }

   // In-place addition

   /*
   * Vector-vector in-place addition, a[i] += b[i] (real).
   */
   void addEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] += b[i];
      }
   }

   /*
   * Vector-scalar in-place addition, a[i] += b (real).
   */
   void addEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] += b;
      }
   }


   // In-place subtraction

   /*
   * Vector-vector in-place subtraction, a[i] -= b[i] (real).
   */
   void subEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] -= b[i];
      }
   }

   /*
   * Vector-scalar in-place subtraction, a[i] -= b (real).
   */
   void subEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] -= b;
      }
   }

   // In-place multiplication

   /*
   * Vector-vector in-place multiplication, a[i] *= b[i] (real).
   */
   void mulEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] *= b[i];
      }
   }

   /*
   * Vector-scalar in-place multiplication, a[i] *= b (real).
   */
   void mulEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] *= b;
      }
   }

   // In-place division

   /*
   * Vector-vector in-place division, a[i] /= b[i] (real).
   */
   void divEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] /= b[i];
      }
   }

   /*
   * Vector-scalar in-place division, a[i] /= b (real).
   */
   void divEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] /= b;
      }
   }

   // Exponentiation

   /*
   * Vector exponentiation, a[i] = exp(b[i]) (real).
   */
   void expV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = exp(b[i]);
      }
   }

   // Square

   /*
   * Vector elementwise square, a[i] = b[i]*b[i] (real).
   */
   void sqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i]*b[i];
      }
   }

   // Absolute magnitude

   /*
   * Vector elementwise square, a[i] = b[i]*b[i] (real).
   */
   void absV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = std::abs(b[i]);
      }
   }

   // Combined operations 

   /*
   * Vector in-place addition with coefficient, a[i] += b[i]*c (real).
   */
   void addEqVc(Array<double>& a, Array<double> const & b, const double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] += b[i]*c;
      }
   }

} // namespace VecOp
} // namespace Pscf

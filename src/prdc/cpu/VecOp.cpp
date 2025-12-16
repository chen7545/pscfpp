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
   * Vector assignment, a[i] = b[i] (real).
   */
   template <>
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
   * Vector assignment to a scalar, a[i] = b (real).
   */
   template <>
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
   template <>
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
   template <>
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
   template <>
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
   template <>
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
   template <>
   void mulVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i]*c[i];
      }
   }

   /*
   * Vector-scalar multiplication, a[i] = b[i] * c (real).
   */
   template <>
   void mulVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i]*c;
      }
   }

   // Division

   /*
   * Vector-vector division, a[i] = b[i] / c[i] (real).
   */
   template <>
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
   template <>
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
   template <>
   void divSV(Array<double>& a, double b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(c.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = b/c[i];
      }
   }

   // Exponentiation

   /*
   * Vector exponentiation, a[i] = exp(b[i]) (real).
   */
   template <>
   void expV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);
      for (int i = 0; i < n; ++i) {
         a[i] = exp(b[i]);
      }
   }

   // Compound addition

   /*
   * Vector-vector in-place addition, a[i] += b[i].
   */
   template <>
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
   * Vector-scalar in-place addition, a[i] += b.
   */
   template <>
   void addEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] += b;
      }
   }


   // Compound subtraction

   /*
   * Vector-vector in-place subtraction, a[i] -= b[i].
   */
   template <>
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
   * Vector-scalar in-place subtraction, a[i] -= b.
   */
   template <>
   void subEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] -= b;
      }
   }

   // Compound multiplication

   /*
   * Vector-vector in-place multiplication, a[i] *= b[i].
   */
   template <>
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
   * Vector-scalar in-place multiplication, a[i] *= b.
   */
   template <>
   void mulEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] *= b;
      }
   }

   // Compound division

   /*
   * Vector-vector in-place division, a[i] /= b[i].
   */
   template <>
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
   * Vector-scalar in-place division, a[i] /= b.
   */
   template <>
   void divEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] /= b;
      }
   }

} // namespace VecOp
} // namespace Pscf


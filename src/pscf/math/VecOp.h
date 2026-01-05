#ifndef PSCF_VEC_OP_H
#define PSCF_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Element-wise vector operations on real-valued Cpu arrays.
   *
   * Operations that are performed by these functions include addition,
   * subtraction, multiplication, division, exponentiation, and assignment.
   * The function names will, correspondingly, begin with "add", "sub",
   * "mul", "div", "exp", or "eq" to indicate the relevant operation.
   * Functions are also included to perform compound assignment operations,
   * i.e.  those that are performed using +=, -=, *=, and /= in C++. These
   * functions have names that begin with "addEq", "subEq", "mulEq", and
   * "divEq", respectively.
   *
   * The output (the LHS of the vector operation) is always the first
   * parameter passed to the function. The input argument(s) (on the RHS
   * of the vector operation) may be vectors or scalars. If an argument 
   * is a vector (scalar), the function name will contain a V (S). For
   * example, the function addVV(A,B,C) implements vector-vector addition
   * A[i] = B[i] + C[i], while addVS(A,B,c) implements vector-scalar
   * addition A[i] = B[i] + c in which c is a scalar that is added to 
   * every element of B. In commutative operations involving both vectors 
   * and scalars, the vectors are listed first. So, for example, addVS 
   * exists, but addSV does not.
   *
   * \ingroup Prdc_Cpu_Module
   * @{
   */
   namespace VecOp {

      // Assignment

      /**
      * Vector assignment, a[i] = b[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void eqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector assignment, a[i] = b (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void eqS(Array<double>& a, double b);

      // Addition

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (real)
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void addVV(Array<double>& a, 
                 Array<double> const & b, 
                 Array<double> const & c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void addVS(Array<double>& a, Array<double> const & b, double c);

      // Subtraction

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (real)
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void subVV(Array<double>& a,
                 Array<double> const & b, Array<double> const & c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void subVS(Array<double>& a, Array<double> const & b, double c);

      // Multiplication

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void mulVV(Array<double>& a,
                 Array<double> const & b, Array<double> const & c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void mulVS(Array<double>& a, Array<double> const & b, double c);

      // Division

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void divVV(Array<double>& a,
                 Array<double> const & b, Array<double> const & c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void divVS(Array<double>& a, Array<double> const & b, double c);

      /**
      * Vector division, a[i] = b / c[i].
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param c  input array (RHS)
      */
      void divSV(Array<double>& a, double b, Array<double> const & c);

      // In-place addition

      /**
      * Vector-vector in-place addition, a[i] += b[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void addEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place addition, a[i] += b (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void addEqS(Array<double>& a, double b);

      // Compound subtraction

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void subEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar subtraction in-place, a[i] -= b (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void subEqS(Array<double>& a, double b);

      // In-place multiplication

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void mulEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void mulEqS(Array<double>& a, double b);

      // Compound division

      /**
      * Vector-vector in-place division, a[i] /= b[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void divEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place division, a[i] /= b.
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void divEqS(Array<double>& a, double b);

      // Exponentiation

      /**
      * Vector exponentiation, a[i] = exp(b[i]) (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void expV(Array<double>& a, Array<double> const & b);

      // Square

      /**
      * Element-wise square, a[i] = b[i]*b[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void sqV(Array<double>& a, Array<double> const & b);

      // Absolute magnitude

      /**
      * Element-wise absolute magnitude, a[i] = abs(b[i]) (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void absV(Array<double>& a, Array<double> const & b);

   /** @} */

   } // namespace VecOp

} // namespace Pscf

#endif

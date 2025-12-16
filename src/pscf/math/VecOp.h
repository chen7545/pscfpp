#ifndef PSCF_MATH_VEC_OP_H
#define PSCF_MATH_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   /**
   * Functions that perform element-wise vector operations on vectors.
   *
   * Operations that are performed by these functions include addition,
   * subtraction, multiplication, division, exponentiation, and assignment.
   * The function names, correspondingly, begin with "add", "sub", "mul",
   * "div", "exp", or "eq" to indicate the relevant operation. Functions 
   * that perform "in-place" arithmetic assignment operations, analogous 
   * to the C/C++ operators +=, -=, *=, and /=, have names that begin 
   * with "addEq", "subEq", "mulEq", and "divEq".
   *
   * The output (the LHS of a vector operation) is always the first
   * parameter passed to the function. The input argument(s) (on the RHS
   * of the vector operation) may be vectors or scalars. If an argument is
   * a vector (or scalar), the function name will contain a V (or S). For
   * example, the function addVV(A,B,C) implements vector-vector addition
   * A[i] = B[i] + C[i], while addVS(A,B,c) implements vector-scalar
   * addition A[i] = B[i] + c in which c is a scalar that is added to every
   *  element of B. In commutative operations involving both vectors and
   * scalars, the vectors are listed first. So, for example, addVS exists,
   * but addSV does not.
   *
   * \ingroup Prdc_Cpu_Module
   * @{
   */
   namespace VecOp {

      // Assignment

      /**
      * Vector-vector assignment, a[i] = b[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void eqV(AT<VT>& a, AT<VT> const & b);

      /**
      * Vector-vector assignment, a[i] = b[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void eqV(AT<VT>& a, AT<RT> const & b);

      /**
      * Vector-scalar assignment, a[i] = b .
      *
      * \param a  array (LHS)
      * \param b  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void eqS(AT<VT>& a, VT b);

      /**
      * Vector-scalar assignment, a[i] = b (mixed).
      *
      * \param a  array (LHS)
      * \param b  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void eqS(AT<VT>& a, RT b);

      // Addition

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void addVV(AT<VT>& a, AT<VT> const & b, AT<VT> const & c);

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void addVV(AT<VT> & a, AT<VT> const & b, AT<RT> const & c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void addVS(AT<VT> & a, AT<VT> const & b, VT c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void addVS(AT<VT> & a, AT<VT> const & b, RT c);

      // Subtraction

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void subVV(AT<VT> & a, AT<VT> const & b, AT<VT> const & c);

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void subVV(AT<VT> & a, AT<VT> const & b, AT<RT> const & c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void subVS(AT<VT> & a, AT<VT> const & b, VT c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void subVS(AT<VT> & a, AT<VT> const & b, RT c);

      // Multiplication

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void mulVV(AT<VT> & a, AT<VT> const & b, AT<VT> const & c);

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void mulVV(AT<VT> & a, AT<VT> const & b, AT<RT> const & c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void mulVS(AT<VT> & a, AT<VT> const & b, VT c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void mulVS(AT<VT> & a, AT<VT> const & b, RT c);

      // Division

      /**
      * Vector-vector division, a[i] = b[i] / c[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void divVV(AT<VT> & a,
                 AT<VT> const & b,
                 AT<VT> const & c);

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void divVV(AT<VT> & a, AT<VT> const & b, AT<RT> const & c);

      /**
      * Vector-scalar division, a[i] = b[i] / c .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void divVS(AT<VT> & a, AT<VT> const & b, VT c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (mixed).
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      * \param c  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void divVS(AT<VT> & a, AT<VT> const & b, RT c);

      /**
      * Vector-scalar division, a[i] = b / c[i] .
      *
      * \param a  array (LHS)
      * \param b  scalar (RHS)
      * \param c  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void divSV(AT<VT> & a, VT b, AT<VT> const & c);

      // In-place addition

      /*
      * Vector-vector in-place addition, a[i] += b[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void addEqV(AT<VT> & a, AT<VT> const & b);

      /*
      * Vector-vector in-place addition, a[i] += b[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void addEqV(AT<VT> & a, AT<RT> const & b);

      /*
      * Vector-scalar in-place addition, a[i] += b .
      *
      * \param a  array (LHS)
      * \param b  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void addEqS(AT<VT> & a, VT b);

      /**
      * Vector-scalar in-place addition, a[i] += b (mixed).
      *
      * \param a  array (LHS)
      * \param b  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void addEqS(AT<VT> & a, RT b);

      // In-place subtraction

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void subEqV(AT<VT>& a, AT<VT> const & b);

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void subEqV(AT<VT>& a, AT<RT> const & b);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b .
      *
      * \param a  array (LHS)
      * \param b  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void subEqS(AT<VT>& a, VT b);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (mixed).
      *
      * \param a  array (LHS)
      * \param b  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void subEqS(AT<VT>& a, RT b);

      // In-place multiplication

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void mulEqV(AT<VT>& a, AT<VT> const & b);

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void mulEqV(AT<VT>& a, AT<RT> const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b[i] .
      *
      * \param a  array (LHS)
      * \param b  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void mulEqV(AT<VT>& a, VT const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b (mixed)
      *
      * \param a  array (LHS)
      * \param b  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void mulEqS(AT<VT>& a, RT b);

      // In-place division

      /**
      * Vector-vector in-place division, a[i] /= b[i] .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void divEqV(AT<VT>& a, AT<VT> const & b);

      /**
      * Vector-vector in-place division, a[i] /= b[i] (mixed).
      *
      * \param a  array (LHS)
      * \param b  real array (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void divEqV(AT<VT>& a, AT<RT> const & b);

      /**
      * Vector-scalar in-place division, a[i] /= b .
      *
      * \param a  array (LHS)
      * \param b  scalar (RHS)
      */
      template <template <typename> class AT, typename VT>
      void divEqS(AT<VT>& a, VT b);

      /**
      * Vector-scalar in-place division, a[i] /= b (mixed).
      *
      * \param a  array (LHS)
      * \param b  real scalar (RHS)
      */
      template <template <typename> class AT, typename VT, typename RT>
      void divEqS(AT<VT>& a, RT b);

      // Exponentiation

      /**
      * Vector exponentiation, a[i] = exp(b[i]) .
      *
      * \param a  array (LHS)
      * \param b  array (RHS)
      */
      template <template <typename> class AT, typename VT>
      void expV(AT<VT> & a, AT<VT> const & b);

   /** @} */

   } // namespace VecOp

} // namespace Pscf

#endif

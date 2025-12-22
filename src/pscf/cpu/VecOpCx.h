#ifndef PSCF_VEC_OP_CX_H
#define PSCF_VEC_OP_CX_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Functions that perform element-wise vector operations on the Cpu.
   *
   * Operations that are performed by these functions include assignment,
   * addition, subtraction, multiplication, division, and exponentiation.
   * The function names will, correspondingly, begin with "eq", "add",
   * "sub", "mul", "div", or "exp" to indicate the relevant operation.
   * Functions that perform arithmetic "in-place" assignment operations,
   * analogous to those performed using +=, -=, *=, and /= operators, 
   * have names that begin with "addEq", "subEq", "mulEq", and "divEq".
   *
   * The output array (the LHS of a vector operation) is always the first
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
      * Vector-vector assignment, a[i] = b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void eqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector-vector assignment, a[i] = b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void eqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar assignment, a[i] = b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void eqS(Array<fftw_complex>& a, fftw_complex b);

      /**
      * Vector-scalar assignment, a[i] = b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void eqS(Array<fftw_complex>& a, double b);

      // Real and imaginary parts

      /**
      * Copy real part of a complex array to a real array.
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      */
      void real(Array<double>& a, Array<fftw_complex> const & b);

      /**
      * Copy imaginary part of a complex array to a real array.
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      */
      void imag(Array<double>& a, Array<fftw_complex> const & b);

      // Addition

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void addVV(Array<fftw_complex>& a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void addVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void addVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void addVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 double c);

      // Subtraction

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void subVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void subVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void subVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 double c);

      // Multiplication

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] (complex)
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void mulVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] (mixed)
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void mulVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (complex)
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void mulVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void mulVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 double c);

      // Division

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void divVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void divVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void divVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void divVS(Array<fftw_complex> & a,
                       Array<fftw_complex> const & b,
                       double c);

      // Exponentiation

      /**
      * Vector exponentiation, a[i] = exp(b[i]) (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void expV(Array<fftw_complex> & a, Array<fftw_complex> const & b);

      // In-place addition

      /*
      * Vector-vector in-place addition, a[i] += b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void addEqV(Array<fftw_complex> & a, Array<fftw_complex> const & b);

      /*
      * Vector-vector in-place addition, a[i] += b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void addEqV(Array<fftw_complex> & a, Array<double> const & b);

      /*
      * Vector-scalar in-place addition, a[i] += b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void addEqS(Array<fftw_complex> & a, fftw_complex b);

      /**
      * Vector-scalar in-place addition, a[i] += b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void addEqS(Array<fftw_complex> & a, double b);

      // In-place subtraction

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void subEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void subEqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void subEqS(Array<fftw_complex>& a, fftw_complex b);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void subEqS(Array<fftw_complex>& a, double b);

      // In-place multiplication

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void mulEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void mulEqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void mulEqV(Array<fftw_complex>& a, fftw_complex const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b (mixed)
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void mulEqS(Array<fftw_complex>& a, double b);

      // In-place division

      /**
      * Vector-vector in-place division, a[i] /= b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void divEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector-vector in-place division, a[i] /= b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void divEqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place division, a[i] /= b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void divEqS(Array<fftw_complex>& a, fftw_complex b);

      /**
      * Vector-scalar in-place division, a[i] /= b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void divEqS(Array<fftw_complex>& a, double b);

   /** @} */

   } // namespace VecOp

} // namespace Pscf

#endif

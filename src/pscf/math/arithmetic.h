#ifndef PSCF_MATH_ARITHMETIC_H
#define PSCF_MATH_ARITHMETIC_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <cmath>
#include <complex>
//#include <iostream>

namespace Pscf {

   /**
   * \defgroup Pscf_Math_Arithmetic_Module Arithmetic Functions
   *
   * Declarations of function templates for arithmetic operations for 
   * real and complex data types. Definitions are given for some explicit 
   * instantiations involving only real data types (double and float), 
   * but only declarations are given in this file for functions that 
   * involve complex types.  Definitions of instantiations that involve
   * involve complex types are given for the complex data types used by 
   * the FFT and cufft FFT libraries in files named complex.h located in 
   * directories src/prdc/cpu and src/prdc/cuda, respectively.
   *
   * Convention: Throughout, a template argument T represents a numerical 
   * type that may be either real or complex, a template argument CT 
   * represents a complex data type, while RT represents a corresponding 
   * real data type.
   *
   * Convention: Functions for which the result or output may be a
   * complex number provide this as a modified value of the first parameter
   * of the function, which must be passed as a non-const reference.
   * Functions for which the output must always be real value (such as an
   * absolute value function) instead provide this as the function return
   * value. This convention allows for the use of C data types for which
   * no assignment operator is defined to represent complex numbers.
   *
   * \ingroup Pscf_Math_Module
   */

   /**
   * Declaration of trait class associated with a complex data type.
   *
   * Instantiations of this class must provide a typedef named
   * Real that is name of the real type used for the real and 
   * imaginary parts.
   */
   template <typename CT> class complexTrait;

   /*
   * The remainder of this file contains declarations of function
   * templates for complex arithmetic that all belong to doxygen topic
   * module Pscf_Math_Arithmetic_Module, which is documented above.
   */

   // Real and imaginary components of complex numbers

   /**
   * Return the real part of a complex number.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex argument (input)
   */
   template <typename CT, typename RT>
   RT real(CT const & a);

   /**
   * Return the imaginary part of a complex number.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex argument (input)
   */
   template <typename CT, typename RT>
   RT imag(CT const & a);

   // Absolute magnitude

   /**
   * Return absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex argument (in)
   */
   template <typename CT, typename RT>
   RT abs(CT const & a);

   /**
   * Return square of absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex argument (in)
   */
   template <typename CT, typename RT>
   RT absSq(CT const & a);

   // Complex Conjugation

   /**
   * Compute complex conjugate, z = a^*.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex conjugate of argument (out)
   * \param a complex argument (in)
   */
   template <typename CT>
   void conj(CT & z, CT const & a);

   /**
   * In place complex conjugation of a complex number, a = a^* .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and complex conjugate (out)
   */
   template <typename CT>
   void conj(CT & a);

   // Assignment

   /**
   * Assign a value from an input of the same type.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   template <typename T>
   void assign(T & z, T const & a);

   /**
   * Assign a double value.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   template <> inline
   void assign(double & z, double const & a)
   {  z = a; }

   /**
   * Assign a real float value.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   template <> inline
   void assign(float & z, float const & a)
   {  z = a; }

   /**
   * Assign one std::complex<RT> variable to another.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   template <typename RT> 
   void assign(std::complex<RT> & z, std::complex<RT> const & a)
   {  z = a; }

   /**
   * Create a complex number from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   template <typename CT, typename RT>
   void assign(CT & z, RT const & a, RT const & b);

   /**
   * Create std::complex<RT> from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   template <typename RT>
   void assign(std::complex<RT> & z, RT const & a, RT const & b)
   {  z = std::complex<RT>(a, b); }

   /**
   * Assign a real input to a complex variable.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   template <typename CT, typename RT>
   void assign(CT & z, RT const & a);

   /**
   * Assign a real input to a std::complex variable.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   template <typename RT>
   void assign(std::complex<RT> & z, RT const & a)
   {  z = a; }

   /**
   * Assign a std::complex input to a complex CT variable, z=a.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a std::complex (in)
   */
   template <typename CT, typename RT>
   void assign(CT & z, std::complex<RT> const & a);

   /**
   * Assign a complex CT input to a std::complex variable, z=a.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z std::complex<RT> (out)
   * \param a complex (in)
   */
   template <typename CT, typename RT>
   void assign(std::complex<RT> & z, CT const & a);

   // Addition

   /**
   * Addition of two numbers of the same type, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   template <typename T>
   void add(T & z, T const & a, T const & b);

   /**
   * Addition of two double precision numbers, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   template <>  inline
   void add(double & z, double const & a, double const & b)
   {  z = a + b; }

   /**
   * Addition of two real float numbers, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   template <> inline
   void add(float & z, float const & a, float const & b)
   {  z = a + b; }

   /**
   * Addition of two std::complex variables, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   template <typename RT> 
   void add(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a + b; }

   /**
   * Addition of a complex and real number, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   template <typename CT, typename RT>
   void add(CT & z, CT const & a, RT const & b);

   /**
   * Addition of a std::complex and real number, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   template <typename RT>
   void add(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z += a; }

   /**
   * In place addition, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <typename T>
   void addEq(T & a, T const & b);

   /**
   * In place addition of double precision real numbers, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <> inline
   void addEq(double & a, double const & b)
   {  a += b; }

   /**
   * In place addition of double precision real numbers, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <> inline
   void addEq(float & a, float const & b)
   {  a += b; }

   /**
   * In place addition of std::complex numbers, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <typename RT> 
   void addEq(std::complex<RT> & a, std::complex<RT> const & b)
   {  a += b; }

   /**
   * In place addition of a complex and real number, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b real summand (in)
   */
   template <typename CT, typename RT>
   void addEq(CT & a, RT const & b);

   /**
   * In place addition of std::complex and real numbers, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <typename RT> 
   void addEq(std::complex<RT> & a, RT const & b)
   {  a += b; }

   // Subtraction

   /**
   * Subtraction, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   template <typename T>
   void sub(T & z, T const & a, T const & b);

   /**
   * Subtraction of double precision numbers, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   template <>  inline
   void sub(double & z, double const & a, double const & b)
   {  z = a - b; }

   /**
   * Subtraction of real float numbers, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   template <> inline
   void sub(float & z, float const & a, float const & b)
   {  z = a - b; }

   /**
   * Subtraction of std::complex numbers, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   template <typename RT> inline
   void sub(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a - b; }

   /**
   * Subtraction of a real number from a complex number, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   template <typename CT, typename RT>
   void sub(CT & z, CT const & a, RT const & b);

   /**
   * Subtraction of a real number from a complex number, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   template <typename RT>
   void sub(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z = a - b; }

   /**
   * In place subtraction of two numbers, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <typename T>
   void subEq(T & a, T const & b);

   /**
   * In place subtraction of double real numbers, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <> inline
   void subEq(double & a, double const & b)
   {  a -= b; }

   /**
   * In place subtraction of float real numbers, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <> inline
   void subEq(float & a, float const & b)
   {  a -= b; }

   /**
   * In place subtraction of std::complex numbers, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <typename RT> 
   void subEq(std::complex<RT> & a, std::complex<RT> const & b)
   {  a -= b; }

   /**
   * In place subtraction of real number from a complex number, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b real argument (in)
   */
   template <typename CT, typename RT>
   void subEq(CT & a, RT const & b);

   /**
   * In place subtraction of real from std::complex number, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <typename RT> 
   void subEq(std::complex<RT> & a, RT const & b)
   {  a -= b; }

   /**
   * Return square of the absolute magnitude of a complex difference.
   *
   * This function returns |a-b|^2 for complex a and b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   template <typename CT, typename RT>
   RT absSqDiff(CT const & a, CT const & b);

   // Multiplication

   /**
   * Multiplication of two numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   template <typename T>
   void mul(T & z, T const & a, T const & b);

   /**
   * Multiplication of double real numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   template <> inline
   void mul(double & z, double const & a, double const & b)
   {  z = a * b; }

   /**
   * Multiplication of float real numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   template <> inline
   void mul(float & z, float const & a, float const & b)
   {  z = a * b; }

   /**
   * Multiplication of std::complex numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   template <typename RT> 
   void mul(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a * b; }

   /**
   * Multiplication of complex and real numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   template <typename CT, typename RT>
   void mul(CT & z, CT const & a, RT const & b);

   /**
   * Multiplication of std::complex and real numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   template <typename RT>
   void mul(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z = a * b; }

   /**
   * In place multiplication, a *= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a factor (in) and product (out)
   * \param b factor (in)
   */
   template <typename T>
   void mulEq(T & a, T const & b);

   /**
   * In place multiplication of double real numbers, a *= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a factor (in) and product (out)
   * \param b factor (in)
   */
   template <> inline
   void mulEq(double & a, double const & b)
   {  a *= b; }

   /**
   * In place multiplication of float real numbers, a *= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a factor (in) and product (out)
   * \param b factor (in)
   */
   template <> inline
   void mulEq(float & a, float const & b)
   {  a *= b; }

   /**
   * In place multiplication of a complex and real number, a *= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b real factor (in)
   */
   template <typename CT, typename RT>
   void mulEq(CT & a, RT const & b);

   /**
   * Compute complex square of a complex number, z = a * a.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   template <typename CT>
   void square(CT & z, CT const & a);

   /**
   * Compute complex square of a std::complex number, z = a * a.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   template <typename RT>
   void square(std::complex<RT> & z, std::complex<RT> const & a)
   {  z = a * a; }

   // Division

   /**
   * Division of two numbers, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   template <typename T>
   void div(T & z, T const & a, T const & b);

   /**
   * Division of two double precision real numbers, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   template <> inline
   void div(double & z, double const & a, double const & b)
   {  z = a/b; }

   /**
   * Division of two float real numbers, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   template <> inline
   void div(float & z, float const & a, float const & b)
   {  z = a / b; }

   /**
   * Division of two std::complex numbers, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   template <typename RT> 
   void div(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a / b; }

   /**
   * Division of a complex number by a real number, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   template <typename CT, typename RT>
   void div(CT & z, CT const & a, RT const & b);

   /**
   * Division of a std::complex number by a real number, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   template <typename RT>
   void div(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z = a / b; }

   /**
   * In place division of two numbers, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a numerator (in) and ratio (out)
   * \param b denominator (in)
   */
   template <typename T>
   void divEq(T & a, T const & b);

   /**
   * In place division of double real numbers, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   template <> inline
   void divEq(double & a, double const & b)
   {  a /= b; }

   /**
   * In place division of float real numbers, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   template <> inline
   void divEq(float & a, float const & b)
   {  a /= b; }

   /**
   * In place division of std::complex numbers, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   template <typename RT> 
   void divEq(std::complex<RT> & a, std::complex<RT> const & b)
   {  a /= b; }

   /**
   * In place division of a complex number by a real number, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   template <typename CT, typename RT>
   void divEq(CT & a, RT const & b);

   /**
   * In place division of std::complex number by a real number, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   template <typename RT>
   void divEq(std::complex<RT> & a, RT const & b)
   {  a /= b; }

   // Inversion

   // Inversion

   /**
   * Inverse of a number, z = 1 / a .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   template <typename T>
   void inverse(T & z, T const & a);

   /**
   * Inverse of a double precision real number, z = 1/a.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   template <> inline
   void inverse(double & z, double const & a)
   {  z = 1.0/a; }

   /**
   * Inverse of a  float real number, z = 1 / a .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   template <> inline
   void inverse(float & z, float const & a)
   {  z = 1.0/a; }

   // Exponential function

   /**
   * Exponentiation of a number, z = exp(a).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   template <typename T>
   void assignExp(T & z, T const & a);

   /**
   * Exponentiation of a real double number, z = exp(a).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   template <> inline
   void assignExp(double & z, double const & a)
   {  z = std::exp(a); }

   /**
   * Exponent of a real float number, z = exp(a).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   template <> inline
   void assignExp(float & z, float const & a)
   {  z = std::exp(a); }

   // Natural logarithm function

   /**
   * Logarithm of a number, z = exp(a) (base template).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   template <typename T>
   void assignLog(T & z, T const & a);

   /**
   * Logarithm of a real double number, z = exp(a).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   template <> inline
   void assignLog(double & z, double const & a)
   {  z = std::log(a); }

   /**
   * Logarithm of a real float number, z = exp(a).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   template <> inline
   void assignLog(float & z, float const & a)
   {  z = std::log(a); }

} // namespace Pscf
#endif

#ifndef PSCF_CUDA_VEC_OP_H
#define PSCF_CUDA_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>

namespace Pscf {

/**
* Functions that perform element-wise vector operations on the GPU.
*
* The functions defined in this file are wrappers for CUDA kernels that
* perform the actual vector operations. The kernels themselves are only
* intended to be called through their wrappers, so the kernel are 
* defined in an anonymous namespace in VecOp.cu.
*
* The functions declared in this header operate on DeviceArray objects 
* with elements of type cudaReal or cudaComplex. The operations that are
* performed by these functions include assignment, addition, subtraction,
* multiplication, division, and exponentiation. Function names,
* correspondingly, begin with "eq", "add", "sub", "mul", "div", or
* "exp" to indicate the operation being performed. Functions that
* perform in-place arithmetic assignment operations, which are 
* analogous to those performed using +=, -=, *=, and /= in C++, have
* names that begin with "addEq", "subEq", "mulEq", and "divEq".
*
* The functions are overloaded to perform their respective operations
* on any combination of cudaReal and cudaComplex input arrays, except
* those that would result in division by a complex number.
*
* The output (the LHS of the vector operation) is always the first
* parameter passed to the function. The input argument(s) (on the RHS
* of the vector operation) may be vectors or scalars. If an argument is
* a vector (scalar), the function name will contain a V (S). For example,
* addVV(A,B,C) implements vector-vector addition A[i] = B[i] + C[i],
* while addVS(A,B,c) implements vector-scalar addition A[i] = B[i] + c
* in which c is a scalar that is added to every element of B. In
* commutative binary operations involving a vector and a scalar, the
* vector is listed first. So, for example, addVS exists, but addSV does
* not.
*
* Two wrapper functions are provided for each vector operation:
* - The first accepts only the output array and the necessary input
*   arrays / scalars. In these functions, each input array must be at
*   least as long as the output array, and the element-wise operation
*   will be performed for every element of the output array. All
*   arrays will be indexed starting at element 0.
* - The second allows for vector operations to be performed using only
*   subsections (slices) of the input and output arrays. These functions
*   require additional parameters: one index for each array involved in
*   the operation (output and input), representing the element of each
*   array at which to begin the slice, and an integer n, representing
*   the size of the slices. Before calling the CUDA kernel, these
*   functions check to ensure that the slices do not contain any indices
*   that exceed the length of the corresponding arrays.
*
* Additional functions that perform multiple operations within a single
* kernel are defined in VecOpMisc. This collection is not comprehensive
* and is added to as-needed during the development of this software.
* VecOpMisc.h is included at the end of VecOp.h so that any code that
* includes VecOp.h will also include VecOpMisc.h.
*
* \ingroup Pscf_Cuda_Module
* @{
*/
namespace VecOp {

// Assignment operations:
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector assignment, a[i] = b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void eqV(DeviceArray<cudaReal>& a, 
         DeviceArray<cudaReal> const & b,
         const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
*/
inline
void eqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  eqV(a, b, 0, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b[i], (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void eqV(DeviceArray<cudaComplex>& a, 
         DeviceArray<cudaComplex> const & b,
         const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], (complex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline
void eqV(DeviceArray<cudaComplex>& a,
                DeviceArray<cudaComplex> const & b)
{  eqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-scalar assignment, a[i] = b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void eqS(DeviceArray<cudaReal>& a, const cudaReal b,
         const int beginIdA, const int n);

/**
* Vector-scalar assignment, a[i] = b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
*/
inline
void eqS(DeviceArray<cudaReal>& a, const cudaReal b)
{  eqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar assignment, a[i] = b, (complex).
*
* \param a  complex array (LHS)
* \param b  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void eqS(DeviceArray<cudaComplex>& a, 
         const cudaComplex b,
         const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, (complex).
*
* \param a  complex array (LHS)
* \param b  complex scalar (RHS)
*/
inline
void eqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
{  eqS(a, b, 0, a.capacity()); }

// Addition operations
// ~~~~~~~~~~~~~~~~~~~

/**
* Vector-vector addition, a[i] = b[i] + c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaReal>& a,
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c,
           const int beginIdA, const int beginIdB, const int beginIdC,
           const int n);

/**
* Vector-vector addition, a[i] = b[i] + c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
*/
inline
void addVV(DeviceArray<cudaReal>& a,
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector addition, a[i] = b[i] + c[i], (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaComplex> const & c,
           const int beginIdA, const int beginIdB, const int beginIdC,
           const int n);

/**
* Vector-vector addition, a[i] = b[i] + c[i], (complex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline
void addVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaComplex> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector addition, a[i] = b[i] + c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaComplex> & a,
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaComplex> const & c,
           const int beginIdA, const int beginIdB, const int beginIdC,
           const int n);

/**
* Vector addition, a[i] = b[i] + c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex array (RHS)
*/
inline
void addVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaComplex> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c,
           const int beginIdA, const int beginIdB, const int beginIdC,
           const int n);

/**
* Vector-vector addition, a[i] = b[i] + c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
*/
inline
void addVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-scalar addition, a[i] = b[i] + c (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           const cudaReal c, 
           const int beginIdA, const int beginIdB, 
           const int n);

/**
* Vector-scalar addition, a[i] = b[i] + c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
*/
inline
void addVS(DeviceArray<cudaReal>& a,
           DeviceArray<cudaReal> const & b,
           const cudaReal c)
{  addVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar addition, a[i] = b[i] + c, (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaComplex>& a, 
  	   DeviceArray<cudaComplex> const & b,
           const cudaComplex c, 
	   const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex scalar (RHS)
*/
inline
void addVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b, const cudaComplex c)
{  addVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaReal> const & b,
           const cudaComplex c, 
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex scalar (RHS)
*/
inline
void addVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b,
           const cudaComplex c)
{  addVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar addition, a[i] = b[i] + c (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           const cudaReal c,
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector-scalar addition, a[i] = b[i] + c (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
*/
inline
void addVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           const cudaReal c)
{  addVS(a, b, c, 0, 0, a.capacity()); }


// Subtraction operations
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector-vector subtraction, a[i] = b[i] - c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c, 
	   const int beginIdA, const int beginIdB, const int beginIdC, 
	   const int n);

/**
* Vector-vector subtraction, a[i] = b[i] - c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
*/
inline
void subVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector subtraction, a[i] = b[i] - c[i], (complex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaComplex> const & c, const int beginIdA,
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector-vector subtraction, a[i] = b[i] - c[i], (complex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline
void subVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaComplex> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector subtraction, a[i] = b[i] - c[i] (mixed).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaComplex> const & c, const int beginIdA,
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector-vector subtraction, a[i] = b[i] - c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex array (RHS)
*/
inline
void subVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaComplex> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector subtraction, a[i] = b[i] - c[i] (mixed, c = real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c, const int beginIdA,
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector-vector subtraction, a[i] = b[i] - c[i] (mixed, c = real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real  array (RHS)
*/
inline
void subVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-scalar subtraction, a[i] = b[i] - c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
           const cudaReal c, const int beginIdA,
           const int beginIdB, const int n);

/**
* Vector-scalar subtraction, a[i] = b[i] - c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
*/
inline
void subVS(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           const cudaReal c)
{  subVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar subtraction, a[i] = b[i] - c, (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b,
           const cudaComplex c, 
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector-scalar subtraction, a[i] = b[i] - c, (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex scalar (RHS)
*/
inline
void subVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b, 
           const cudaComplex c)
{  subVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaReal> const & b,
           const cudaComplex c, 
	   const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex scalar (RHS)
*/
inline
void subVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b,
           const cudaComplex c)
{  subVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c (mixed, c = real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b,
           const cudaReal c, 
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c (mixed, c = real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
*/
inline
void subVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           const cudaReal c)
{  subVS(a, b, c, 0, 0, a.capacity()); }


// Multiplication operations
// ~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector-vector multiplication, a[i] = b[i] * c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c, 
	   const int beginIdA, const int beginIdB, const int beginIdC, 
           const int n);

/**
* Vector-vector multiplication, a[i] = b[i] * c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
*/
inline
void mulVV(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector multiplication, a[i] = b[i] * c[i], (complex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaComplex> const & c, 
           const int beginIdA, const int beginIdB, const int beginIdC, 
           const int n);

/**
* Vector-vector multiplication, a[i] = b[i] * c[i], (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex array (RHS)
*/
inline
void mulVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaComplex> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector multiplication, a[i] = b[i] * c[i] (mixed).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
           DeviceArray<cudaComplex> const & c, const int beginIdA,
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector-vector multiplication, a[i]=b[i]*c[i] (mixed).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline
void mulVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaComplex> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector multiplication, a[i]=b[i]*c[i] (mixed, c = real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c, 
           const int beginIdA, const int beginIdB, const int beginIdC, 
           const int n);

/**
* Vector-vector multiplication, a[i]=b[i]*c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
*/
inline
void mulVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-scalar multiplication, a[i] = b[i] * c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaReal>& a,
           DeviceArray<cudaReal> const & b,
           const cudaReal c,
           const int beginIdA, const int beginIdB,
           const int n);

/**
* Vector-scalar multiplication, a[i] = b[i] * c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
*/
inline
void mulVS(DeviceArray<cudaReal>& a,
           DeviceArray<cudaReal> const & b,
           const cudaReal c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar multiplication, a[i] = b[i] * c, (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           const cudaComplex c, 
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector-scalar multiplication, a[i] = b[i] * c, (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  complex scalar (RHS)
*/
inline
void mulVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b, 
	   const cudaComplex c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar multiplication, a[i] = b[i] * c (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param c  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b,
           const cudaComplex c,
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector-scalar multiplication, a[i] = b[i] * c (mixed).
*
* \param a  output array (LHS)
* \param b  real array (RHS)
* \param c  complex scalar (RHS)
*/
inline
void mulVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaReal> const & b, 
	   const cudaComplex c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar multiplication, a[i] = b[i] * c (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           const cudaReal c,
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector-scalar multiplication, a[i] = b[i] * c (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
*/
inline
void mulVS(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaComplex> const & b,
                  const cudaReal c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }


// Division operations
// ~~~~~~~~~~~~~~~~~~~

/**
* Vector-vector division, a[i] = b[i] / c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void divVV(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           DeviceArray<cudaReal> const & c, 
           const int beginIdA, const int beginIdB, const int beginIdC, 
           const int n);

/**
* Vector-vector division, a[i] = b[i] / c[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real array (RHS)
*/
inline
void divVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
                  DeviceArray<cudaReal> const & c)
{  divVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-vector division, a[i] = b[i] / c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void divVV(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c, 
           const int beginIdA, const int beginIdB, const int beginIdC, 
           const int n);

/**
* Vector-vector division, a[i] = b[i] / c[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real array (RHS)
*/
inline
void divVV(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b,
           DeviceArray<cudaReal> const & c)
{  divVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector-scalar division, a[i] = b[i] / c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
           const cudaReal c, const int beginIdA,
           const int beginIdB, const int n);

/**
* Vector-scalar division, a[i] = b[i] / c, (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param c  real scalar (RHS)
*/
inline
void divVS(DeviceArray<cudaReal>& a, 
           DeviceArray<cudaReal> const & b,
           const cudaReal c)
{  divVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector-scalar division, a[i] = b[i] / c (mixed, c real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divVS(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b,
           const cudaReal c, 
           const int beginIdA, const int beginIdB, const int n);

/**
* Vector-scalar division, a[i] = b[i] / c (mixed, c = real).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param c  real scalar (RHS)
*/
inline
void divVS(DeviceArray<cudaComplex>& a,
           DeviceArray<cudaComplex> const & b, 
           const cudaReal c)
{  divVS(a, b, c, 0, 0, a.capacity()); }

/**
* Scalar-vector division, a[i] = b / c[i], (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param c  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void divSV(DeviceArray<cudaReal>& a, 
           const cudaReal b,
           DeviceArray<cudaReal> const & c,
           const int beginIdA, const int beginIdC, const int n);

/**
* Scalar-vector division, a[i] = b / c[i], (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param c  real array (RHS)
*/
inline
void divSV(DeviceArray<cudaReal>& a, 
           const cudaReal b,
           DeviceArray<cudaReal> const & c)
{  divSV(a, b, c, 0, 0, a.capacity()); }

// Exponentiation operations:
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector exponentiation, a[i] = exp(b[i]), (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void expV(DeviceArray<cudaReal>& a,
          DeviceArray<cudaReal> const & b,
          const int beginIdA, const int beginIdB,
          const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
*/
inline
void expV(DeviceArray<cudaReal>& a, 
          DeviceArray<cudaReal> const & b)
{  expV(a, b, 0, 0, a.capacity()); }

/**
* Vector exponentiation, a[i] = exp(b[i]), (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void expV(DeviceArray<cudaComplex>& a,
          DeviceArray<cudaComplex> const & b,
          const int beginIdA, const int beginIdB,
          const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
*/
inline
void expV(DeviceArray<cudaComplex>& a,
          DeviceArray<cudaComplex> const & b)
{  expV(a, b, 0, 0, a.capacity()); }


// Compound operations: addition
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector in-place addition, a[i] += b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place addition, a[i] += b[i] (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
*/
inline
void addEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-vector in-place addition, a[i] += b[i] (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addEqV(DeviceArray<cudaComplex>& a, 
            DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place addition, a[i] += b[i] (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
*/
inline
void addEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaComplex> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-vector in-place addition, a[i] += b[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addEqV(DeviceArray<cudaComplex>& a, 
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place addition, a[i] += b[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
*/
inline
void addEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaReal> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-scalar in-place addition, a[i] += b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void addEqS(DeviceArray<cudaReal>& a, 
            const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place addition, a[i] += b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
*/
inline
void addEqS(DeviceArray<cudaReal>& a, const cudaReal b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar in-place addition, a[i] += b, (complex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void addEqS(DeviceArray<cudaComplex>& a, 
            const cudaComplex b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place addition, a[i] += b, (complex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline
void addEqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar in-place addition, a[i] += b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void addEqS(DeviceArray<cudaComplex>& a, 
            const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place addition, a[i] += b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
*/
inline
void addEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
{  addEqS(a, b, 0, a.capacity()); }


// Compound operations: subtraction
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector in-place subtraction, a[i] -= b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subEqV(DeviceArray<cudaReal>& a, 
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector in-place subtraction, a[i] -= b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
*/
inline
void subEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-vector in-place subtraction, a[i] -= b[i], (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subEqV(DeviceArray<cudaComplex>& a, 
            DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector in-place subtraction, a[i] -= b[i], (complex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline
void subEqV(DeviceArray<cudaComplex>& a,
                   DeviceArray<cudaComplex> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector in-place subtraction, a[i]-=b[i] (mixed).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector in-place subtraction, a[i]-=b[i] (mixed).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline
void subEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaReal> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector in-place subtraction, a[i] -= b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void subEqS(DeviceArray<cudaReal>& a, const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector in-place subtraction, a[i] -= b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
*/
inline
void subEqS(DeviceArray<cudaReal>& a, const cudaReal b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector in-place subtraction, a[i] -= b, (complex).
*
* \param a  complex array (LHS)
* \param b  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void subEqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
            const int beginIdA, const int n);

/**
* Vector in-place subtraction, a[i] -= b, (complex).
*
* \param a  complex array (LHS)
* \param b  complex scalar (RHS)
*/
inline
void subEqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar in-place subtraction, a[i] -= b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void subEqS(DeviceArray<cudaComplex>& a, 
            const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place subtraction, a[i] -= b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
*/
inline
void subEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
{  subEqS(a, b, 0, a.capacity()); }


// Compound operations: multiplication
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector-vector in-place multiplication, a[i] *= b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place multiplication, a[i] *= b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
*/
inline
void mulEqV(DeviceArray<cudaReal>& a,
            DeviceArray<cudaReal> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-vector in-place multiplication, a[i] *= b[i], (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulEqV(DeviceArray<cudaComplex>& a, 
            DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place multiplication, a[i] *= b[i], (complex).
*
* \param a  complex array (LHS)
* \param b  complex array (RHS)
*/
inline
void mulEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaComplex> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-vector in-place multiplication, a[i]*=b[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulEqV(DeviceArray<cudaComplex>& a, 
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place multiplication, a[i] *= b[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
*/
inline
void mulEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaReal> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-scalar in-place multiplication, a[i] *= b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void mulEqS(DeviceArray<cudaReal>& a, const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place multiplication, a[i] *= b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
*/
inline
void mulEqS(DeviceArray<cudaReal>& a, const cudaReal b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar in-place multiplication, a[i] *= b, (complex).
*
* \param a  complex array (LHS)
* \param b  complex scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void mulEqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place multiplication, a[i] *= b, (complex).
*
* \param a  complex array (LHS)
* \param b  complex scalar (RHS)
*/
inline
void mulEqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar in-place multiplication, a[i]*=b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void mulEqS(DeviceArray<cudaComplex>& a, 
            const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place multiplication, a[i]*=b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
*/
inline
void mulEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
{  mulEqS(a, b, 0, a.capacity()); }


// Compound operations: division
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector-vector in-place division, a[i] /= b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divEqV(DeviceArray<cudaReal>& a, 
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place division, a[i] /= b[i], (real).
*
* \param a  real array (LHS)
* \param b  real array (RHS)
*/
inline
void divEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  divEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-vector in-place division, a[i] /= b[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector-vector in-place division, a[i] /= b[i] (mixed).
*
* \param a  complex array (LHS)
* \param b  real array (RHS)
*/
inline
void divEqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaReal> const & b)
{  divEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector-scalar in-place division, a[i] /= b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void divEqS(DeviceArray<cudaReal>& a, 
            const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place division, a[i] /= b, (real).
*
* \param a  real array (LHS)
* \param b  real scalar (RHS)
*/
inline
void divEqS(DeviceArray<cudaReal>& a, const cudaReal b)
{  divEqS(a, b, 0, a.capacity()); }

/**
* Vector-scalar in-place division, a[i] /= b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void divEqS(DeviceArray<cudaComplex>& a, 
            const cudaReal b,
            const int beginIdA, const int n);

/**
* Vector-scalar in-place division, a[i] /= b (mixed).
*
* \param a  complex array (LHS)
* \param b  real scalar (RHS)
*/
inline
void divEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
{  divEqS(a, b, 0, a.capacity()); }

/** @} */

} // namespace VecOp
} // namespace Pscf

// Ensure that if VecOp.h is included, so is VecOpMisc.h
#include "VecOpMisc.h" 
#endif

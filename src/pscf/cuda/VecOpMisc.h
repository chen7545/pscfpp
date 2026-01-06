#ifndef PSCF_CUDA_VEC_OP_MISC_H
#define PSCF_CUDA_VEC_OP_MISC_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace VecOp {

   /*
   * Miscellaneous element-wise vector operations performed on the GPU.
   *
   * Note: this file is included at the end of VecOp.h, so any file that
   * includes VecOp.h will also include this file.
   *
   * The functions defined in this file combine 2 or more element-wise
   * vector operations into a single kernel launch, which will perform
   * the operation faster than consecutively calling multiple of the
   * functions in VecOp.h. This collection of functions is not intended
   * to be comprehensive.  Rather, they are written and included as
   * needed during the development of other code.
   *
   * The names of these functions follow the same conventions as those
   * in VecOp, using add, sub, mul, div, exp, eq, and combinations
   * thereof to indicate the operation(s) being performed. V denotes a
   * a vector, S denotes a scalar, and Vc denotes a vector that is
   * multiplied by a scalar coefficient and then used in another
   * operation. For example, addEqVc(a, b, c) performs a[i] += b[i] * c
   * for all i.
   *
   * Another set of functions defined in this file contain the word Pair,
   * indicating that these functions perform the same operation for a
   * pair of real arrays. For example, eqVPair performs a1[i] = s[i]
   * and a2[i] = s[i] for all i. Performing these operations in pairs
   * is faster because the shared array s only needs to be loaded from 
   * global memory once.
   *
   * A third set of functions defined in this file contain the word
   * "Many", indicating that an undefined number of vectors (>2) are
   * involved in an operation. For example, addVMany adds >2 vectors
   * together by passing an array of vectors, rather than a discrete
   * set of vectors.
   *
   * The functions declared in this file are wrappers for CUDA kernels
   * that perform the actual vector operations. The underlying kernels 
   * are only intended to be called through their wrappers, and so are
   * defined in an anonymous namespace in the file VecOpMisc.cu that
   * is inaccessible outside that file. 
   */

   // Functions that combine multiple VecOp operations

   /**
   * Vector addition w/ coefficient, a[i] = (b[i]*c) + (d[i]*e).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real array 1 (RHS)
   * \param c  input scalar (RHS)
   * \param d  real array 2 (RHS)
   * \param e  input scalar 2 (RHS)
   */
   void addVcVc(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b, cudaReal const c,
                DeviceArray<cudaReal> const & d, cudaReal const e);

   /**
   * 3-vec addition w coeff, a[i] = (b[i]*c) + (d[i]*e) + (f[i]*g).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real array 1 (RHS)
   * \param c  input scalar 1 (RHS)
   * \param d  real array 2 (RHS)
   * \param e  input scalar 2 (RHS)
   * \param f  real array 3 (RHS)
   * \param g  input scalar 3 (RHS)
   */
   void addVcVcVc(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b, cudaReal const c,
                  DeviceArray<cudaReal> const & d, cudaReal const e,
                  DeviceArray<cudaReal> const & f, cudaReal const g);

   /**
   * Vector addition in-place w/ coefficient, a[i] += b[i] * c.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real array (RHS)
   * \param c  input scalar
   */
   void addEqVc(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b,
                cudaReal const c);

   /**
   * Vector subtraction, a[i] = b[i] - c[i] - d.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real array (RHS)
   * \param c  real array (RHS)
   * \param d  input scalar (RHS)
   */
   void subVVS(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               DeviceArray<cudaReal> const & c, cudaReal const d);

   /**
   * Vector division in-place w/ coeff., a[i] /= (b[i] * c).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  complex array (LHS)
   * \param b  real array (RHS)
   * \param c  input scalar (RHS)
   */
   void divEqVc(DeviceArray<cudaComplex>& a,
                DeviceArray<cudaReal> const & b,
                cudaReal const c);

   /**
   * Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real array (RHS)
   * \param c  input scalar
   */
   void expVc(DeviceArray<cudaReal>& a, 
              DeviceArray<cudaReal> const & b,
              cudaReal const c);


   // Pair functions

   /**
   * Vector assignment in pairs, ax[i] = b[i], x = 1, 2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a1  real array 1 (LHS)
   * \param a2  real array 2 (LHS)
   * \param s  shared real array to be assigned to both a1 and a2
   */
   void eqVPair(DeviceArray<cudaReal>& a1, 
                DeviceArray<cudaReal>& a2,
                DeviceArray<cudaReal> const & s);

   /**
   * Vector multiplication in pairs, ax[i] = bx[i] * s[i], x=1,2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a1  real array 1 (LHS)
   * \param a2  real array 2 (LHS)
   * \param b1  real array 1 (RHS)
   * \param b2  real array 2 (RHS)
   * \param s  shared real array to be multiplied by both b1 and b2
   */
   void mulVVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2,
                  DeviceArray<cudaReal> const & b1,
                  DeviceArray<cudaReal> const & b2,
                  DeviceArray<cudaReal> const & s);

   /**
   * In-place vector multiplication in pairs, ax[i] *= s[i], x=1,2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a1  real array 1 (LHS)
   * \param a2  real array 2 (LHS)
   * \param s  shared real array to multiply both a1 and a2 (RHS)
   */
   void mulEqVPair(DeviceArray<cudaReal>& a1, 
                   DeviceArray<cudaReal>& a2,
                   DeviceArray<cudaReal> const & s);


   // Functions of "many" vectors

   /**
   * Add an arbitrary number of vectors pointwise (real).
   *
   * The input array 'vecs' contains the arrays that will be added.
   * The size of array vecs determines the number of vectors that
   * will be added together.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of DeviceArrays to be added (RHS)
   */
   void addVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> > const & vecs);

   /**
   * Add an arbitrary number of vectors pointwise (real).
   *
   * The input array 'vecs' contains const pointers to each array that
   * will be added together. The size of vecs determines the number of
   * vectors that will be added together.
   *
   * This version of addVMany is provided for cases in which one needs to
   * add many arrays that are not already stored together in a DArray. The
   * caller must simply assemble an array of pointers to all of the arrays
   * that should be added, and then pass it to this method.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of pointers to DeviceArrays to be added
   */
   void addVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> const *> const & vecs);

   /**
   * Multiply an undefined number of vectors pointwise (real).
   *
   * The input array 'vecs' contains the arrays that will be multiplied.
   * The size of vecs determines the number of vectors that will be
   * multiplied together.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of DeviceArrays to be multiplied (RHS)
   */
   void mulVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> > const & vecs);

   /**
   * Multiply an undefined number of vectors pointwise.
   *
   * The input array 'vecs' contains const pointers to each array that
   * will be multiplied together. The size of vecs determines the number
   * of vectors that will ultimately be multiplied together pointwise.
   *
   * This version of mulVMany is provided for cases in which one needs to
   * multiply many arrays that are not already stored together in a DArray.
   * The caller must simply assemble an array of pointers to all of the
   * arrays that should be multiplied, and then pass it to this method.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of real arrays to be multiplied (RHS)
   */
   void mulVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> const *> const & vecs);


   // Other useful functions

   /**
   * Squared norm of complex number, a[i] = norm(b[i])^2.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  complex array (RHS)
   */
   void sqNormV(DeviceArray<cudaReal>& a,
                DeviceArray<cudaComplex> const & b);

   /**
   * Norm of complex number to the 4th power, a[i] = norm(b[i])^4.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  complex array (RHS)
   */
   void sqSqNormV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaComplex> const & b);

} // namespace VecOp
} // namespace Pscf
#endif

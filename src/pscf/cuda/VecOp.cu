/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOp.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/cudaErrorCheck.h>
#include <cmath>

namespace Pscf {
namespace VecOp {

// CUDA kernels:
// (defined in anonymous namespace, used within this file only)

namespace {

      /*
      * Vector assignment, a[i] = b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i];
         }
      }

      /*
      * Vector assignment, a[i] = b[i](complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x;
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector assignment, a[i] = b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b;
         }
      }

      /*
      * Vector assignment, a[i] = b(complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqS(cudaComplex* a, const cudaComplex b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b.x;
            a[i].y = b.y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaReal* a,
                  cudaReal const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] + c[i];
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i](complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaComplex* a,
                  cudaComplex const * b,
                  cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c[i].x;
            a[i].y = b[i].y + c[i].y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, b = real).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaComplex* a,
                  cudaReal const * b,
                  cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] + c[i].x;
            a[i].y = c[i].y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, c = real).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaComplex* a,
                  cudaComplex const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c[i];
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaReal* a,
                  cudaReal const * b,
                  const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] + c;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c(complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaComplex* a,
                  cudaComplex const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c.x;
            a[i].y = b[i].y + c.y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c, GPU kernel (mixed, b = real).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaComplex* a,
                  cudaReal const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] + c.x;
            a[i].y = c.y;
         }
      }

      /*
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaComplex* a,
                  cudaComplex const * b,
                  const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c;
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector-vector subtraction, a[i] = b[i] - c[i], GPU kernel (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subVV(cudaReal* a,
                  cudaReal const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] - c[i];
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i](complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVV(cudaComplex* a, cudaComplex const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c[i].x;
            a[i].y = b[i].y - c[i].y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, b = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVV(cudaComplex* a, cudaReal const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] - c[i].x;
            a[i].y = 0.0 - c[i].y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, c = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVV(cudaComplex* a, cudaComplex const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c[i];
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVS(cudaReal* a, cudaReal const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] - c;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c(complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVS(cudaComplex* a, cudaComplex const * b,
                             const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c.x;
            a[i].y = b[i].y - c.y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, b = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVS(cudaComplex* a, cudaReal const * b,
                             const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] - c.x;
            a[i].y = 0.0 - c.y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, c = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVS(cudaComplex* a, cudaComplex const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c;
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVV(cudaReal* a, cudaReal const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] * c[i];
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i](complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVV(cudaComplex* a, cudaComplex const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = (b[i].x * c[i].x) - (b[i].y * c[i].y);
            a[i].y = (b[i].x * c[i].y) + (b[i].y * c[i].x);
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, b = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVV(cudaComplex* a, cudaReal const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] * c[i].x;
            a[i].y = b[i] * c[i].y;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, c = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVV(cudaComplex* a, cudaComplex const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x * c[i];
            a[i].y = b[i].y * c[i];
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVS(cudaReal* a, cudaReal const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] * c;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c(complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVS(cudaComplex* a, cudaComplex const * b,
                             const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = (b[i].x * c.x) - (b[i].y * c.y);
            a[i].y = (b[i].x * c.y) + (b[i].y * c.x);
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, b = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVS(cudaComplex* a, cudaReal const * b,
                             const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] * c.x;
            a[i].y = b[i] * c.y;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, c = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVS(cudaComplex* a, cudaComplex const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x * c;
            a[i].y = b[i].y * c;
         }
      }

      /*
      * Vector division, a[i] = b[i] / c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVV(cudaReal* a, cudaReal const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] / c[i];
         }
      }

      /*
      * Vector division, a[i] = b[i] / c[i], GPU kernel (mixed, c = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVV(cudaComplex* a, cudaComplex const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x / c[i];
            a[i].y = b[i].y / c[i];
         }
      }

      /*
      * Vector division, a[i] = b[i] / c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVS(cudaReal* a, cudaReal const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] / c;
         }
      }

      /*
      * Vector division, a[i] = b[i] / c, GPU kernel (mixed, c = real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVS(cudaComplex* a, cudaComplex const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x / c;
            a[i].y = b[i].y / c;
         }
      }

      /*
      * Vector division, a[i] = b / c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _divSV(cudaReal* a, const cudaReal b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b / c[i];
         }
      }

      // In-place addition

      /*
      * Vector-vector in-place addition, a[i] += b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b[i];
         }
      }

      /*
      * Vector-vector in-place addition, a[i] += b[i] (complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _addEqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i].x;
            a[i].y += b[i].y;
         }
      }

      /*
      * Vector addition in-place, a[i] += b[i], GPU kernel (mixed).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _addEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i];
         }
      }

      /*
      * Vector addition in-place, a[i] += b (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b;
         }
      }

      /*
      * Vector-scalar in-place addition, a[i] += b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqS(cudaComplex* a, const cudaComplex b,
                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b.x;
            a[i].y += b.y;
         }
      }

      /*
      * Vector addition in-place, a[i] += b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b;
         }
      }

      /*
      * Vector subtraction in-place, a[i] -= b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] -= b[i];
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b[i].x;
            a[i].y -= b[i].y;
         }
      }

      /*
      * Vector subtraction in-place, a[i] -= b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b[i];
         }
      }

      /*
      * Vector subtraction in-place, a[i] -= b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] -= b;
         }
      }

      /*
      * Vector subtraction in-place, a[i] -= b, GPU kernel (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqS(cudaComplex* a, const cudaComplex b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b.x;
            a[i].y -= b.y;
         }
      }

      /*
      * Vector subtraction in-place, a[i] -= b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b;
         }
      }

      /*
      * Vector multiplication in-place, a[i] *= b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] *= b[i];
         }
      }

      /*
      * Vector multiplication in-place, a[i] *= b[i], (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaComplex c;
         for (int i = startID; i < n; i += nThreads) {
            c.x = (a[i].x * b[i].x) - (a[i].y * b[i].y);
            c.y = (a[i].x * b[i].y) + (a[i].y * b[i].x);
            a[i].x = c.x;
            a[i].y = c.y;
         }
      }

      /*
      * Vector-vector in-place multiplication, a[i]*=b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x *= b[i];
            a[i].y *= b[i];
         }
      }

      /*
      * Vector-scalar in-place multiplication, a[i] *= b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] *= b;
         }
      }

      /*
      * Vector-scalar multiplication in-place, a[i] *= b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqS(cudaComplex* a, const cudaComplex b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaComplex c;
         for (int i = startID; i < n; i += nThreads) {
            c.x = (a[i].x * b.x) - (a[i].y * b.y);
            c.y = (a[i].x * b.y) + (a[i].y * b.x);
            a[i].x = c.x;
            a[i].y = c.y;
         }
      }

      /*
      * Vector multiplication in-place, a[i] *= b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x *= b;
            a[i].y *= b;
         }
      }

      /*
      * Vector division in-place, a[i] /= b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] /= b[i];
         }
      }

      /*
      * Vector division in-place, a[i] /= b[i], GPU kernel (mixed, b = real).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x /= b[i];
            a[i].y /= b[i];
         }
      }

      /*
      * Vector division in-place, a[i] /= b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] /= b;
         }
      }

      /*
      * Vector division in-place, a[i] /= b, GPU kernel (mixed, b = real).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x /= b;
            a[i].y /= b;
         }
      }

      // Exponentiation

      /*
      * Vector exponentiation, a[i] = exp(b[i]) (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _expV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = exp(b[i]);
         }
      }

      /*
      * Vector exponentiation, a[i] = exp(b[i]) (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _expV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = exp(b[i].x) * cos(b[i].y);
            a[i].y = exp(b[i].x) * sin(b[i].y);
         }
      }

      // Absolute magnitude

      /*
      * Vector absolute magnitude, a[i] = abs(b[i]) (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _absV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = std::fabs(b[i]);
         }
      }

      /*
      * Vector absolute magnitude, a[i] = |b[i]|^2 (complex).
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _absSqV(cudaReal* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal bx, by;
         for (int i = startID; i < n; i += nThreads) {
            bx = b[i].x;
            by = b[i].y;
            a[i] = bx*bx + by*by;
         }
      }

   } // end anonymous namespace

   // CUDA kernel wrappers:

   // Vector assignment, a[i] = b[i] (real).
   void eqV(DeviceArray<cudaReal>& a,
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                  b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector assignment, a[i] = b[i] (complex).
   void eqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                  b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar assignment, a[i] = b (real).
   void eqS(DeviceArray<cudaReal>& a,
            const cudaReal b,
            const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector assignment, a[i] = b (complex).
   void eqS(DeviceArray<cudaComplex>& a,
            const cudaComplex b,
            const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition, a[i] = b[i] + c[i] (real).
   void addVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector addition, a[i] = b[i] + c[i] (complex).
   void addVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition, a[i] = b[i] + c[i] (mixed, b = real).
   void addVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector addition, a[i] = b[i] + c[i] (mixed).
   void addVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar addition, a[i] = b[i] + c (real).
   void addVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar addition, a[i] = b[i] + c (complex).
   void addVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition, a[i] = b[i] + c (mixed, b = real).
   void addVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition, a[i] = b[i] + c (mixed, c = real).
   void addVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c[i] (real).
   void subVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA,
              const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c[i] (cudaComplex).
   void subVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector subtraction, a[i]=b[i]-c[i] (mixed).
   void subVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector subtraction, a[i]=b[i]-c[i] (mixed).
   void subVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c (cudaReal).
   void subVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
               const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c (complex).
   void subVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaComplex c, const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar subtraction, a[i] = b[i] - c (mixed).
   void subVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              const cudaComplex c, const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar subtraction, a[i] = b[i] - c (mixed).
   void subVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i] = b[i] * c[i] (real).
   void mulVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i] = b[i] * c[i] (cudaComplex).
   void mulVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector multiplication, a[i]=b[i]*c[i] (mixed).
   void mulVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector multiplication, a[i]=b[i]*c[i] (mixed).
   void mulVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar multiplication, a[i] = b[i] * c (real).
   void mulVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar multiplication, a[i] = b[i] * c (complex).
   void mulVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar multiplication, a[i] = b[i] * c (mixed).
   void mulVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i] = b[i] * c (mixed).
   void mulVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector division, a[i] = b[i] / c[i] (real).
   void divVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c, const int beginIdA,
              const int beginIdB, const int beginIdC, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector division, a[i] = b[i] / c[i] (mixed).
   void divVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA,
              const int beginIdB, const int beginIdC, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar division, a[i] = b[i] / c (real).
   void divVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector division, a[i] = b[i] / c (mixed, c = real).
   void divVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c, const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Division of scalar by  vector, a[i] = b / c[i] (cudaReal).
   void divSV(DeviceArray<cudaReal>& a,
              const cudaReal b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdC, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divSV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b,
                                    c.cArray() + beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition in-place, a[i] += b[i] (cudaReal).
   void addEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition in-place, a[i] += b[i] (cudaComplex).
   void addEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition in-place, a[i] += b[i] (mixed).
   void addEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar in-place addition, a[i] += b (real).
   void addEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition in-place, a[i] += b (complex).
   void addEqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar in-place addition, a[i] += b (mixed).
   void addEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction in-place, a[i] -= b[i] (cudaReal).
   void subEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector in-place subtraction, a[i] -= b[i] (complex).
   void subEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector in-place subtraction, a[i] -= b[i] (mixed).
   void subEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar in-place subtraction, a[i] -= b (real).
   void subEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar in-place subtraction, a[i] -= b (complex).
   void subEqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar subtraction in-place, a[i] -= b (mixed).
   void subEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // In-place multiplication

   // Vector-vector in-place multiplication, a[i] *= b[i] (real).
   void mulEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication in-place, a[i] *= b[i] (cudaComplex).
   void mulEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication in-place, a[i] *= b[i] (mixed).
   void mulEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication in-place, a[i] *= b (cudaReal).
   void mulEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication in-place, a[i] *= b (cudaComplex).
   void mulEqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication in-place, a[i] *= b (mixed).
   void mulEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // In-place division

   // Vector-vector division in-place, a[i] /= b[i] (real).
   void divEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-vector division in-place, a[i] /= b[i] (mixed).
   void divEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector division in-place, a[i] /= b (real).
   void divEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector division in-place, a[i] /= b (mixed).
   void divEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Exponentiation

   // Vector exponentiation, a[i] = exp(b[i]) (cudaReal).
   void expV(DeviceArray<cudaReal>& a,
             DeviceArray<cudaReal> const & b,
             const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _expV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                   b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector exponentiation, a[i] = exp(b[i]) (cudaComplex).
   void expV(DeviceArray<cudaComplex>& a,
             DeviceArray<cudaComplex> const & b,
             const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _expV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                   b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Absolute magnitudes

   // Vector absolute magnitude, a[i] = abs(b[i]) (real).
   void absV(DeviceArray<cudaReal>& a,
             DeviceArray<cudaReal> const & b,
             const int beginIdA, const int beginIdB,
             const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _absV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                   b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector absolute magnitude squared, a[i] = |b[i]|^2 (complex).
   void absSqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, 
               const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _absSqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

} // namespace VecOp
} // namespace Pscf

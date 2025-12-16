#ifndef PSCF_CUDA_TYPES_H
#define PSCF_CUDA_TYPES_H

#include <cufft.h>

// Toggle single / double precision:
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//#define SINGLE_PRECISION
#define DOUBLE_PRECISION

namespace Pscf {

   /**
   * Complex number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   typedef cufftComplex cudaComplex;
   #else
   #ifdef DOUBLE_PRECISION
   typedef cufftDoubleComplex cudaComplex;
   #endif
   #endif

   /**
   * Real number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   typedef cufftReal cudaReal;
   #else
   #ifdef DOUBLE_PRECISION
   typedef cufftDoubleReal cudaReal;
   #endif
   #endif

}
#endif

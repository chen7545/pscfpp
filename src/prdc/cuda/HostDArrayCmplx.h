#ifndef PRDC_HOST_D_ARRAY_CMPLX_H
#define PRDC_HOST_D_ARRAY_CMPLX_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/HostDArray.h>
#include <prdc/cuda/types.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * HostDArray with complex values. 
   *
   * \ingroup Prdc_Cuda_Module
   */
   class HostDArrayCmplx : public HostDArray<cudaComplex>
   {
   public:

      // Type aliases

      /**
      * Direct base class.
      */
      using Base = HostDArray<cudaComplex>;

      // Type of an element value.
      using Base::ValueType;

      /**
      * Type of real and imaginary parts of an element value.
      */
      using RealType = cudaReal;

      // Public member functions
     
      /**
      * Default constructor.
      */ 
      HostDArrayCmplx();

      /**
      * Allocating constructor.
      *
      * \param n  number of elements for which to allocate memory
      */ 
      HostDArrayCmplx(int n);

      /**
      * Destructor.
      */ 
      ~HostDArrayCmplx();

      // Inherited functions
      using Base::operator=;
      using Base::copySlice;

   };

} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif

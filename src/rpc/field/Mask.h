#ifndef RPC_MASK_H
#define RPC_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Mask.h>       // base class template
#include "FieldIo.h"             // base class template argument
#include <prdc/cpu/RField.h>     // base class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A field to which the total density is constrained.
   *
   * Please refer to the documentation of the base class template Rp::Mask
   * for more complete API documentation for this class template. The
   * public interface of Rpc::Mask is identical to that of the base class
   * template Rp::Mask.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Mask : public Rp::Mask<D, RField<D>, FieldIo<D> >
   {

   public:

      /// Base class typedef
      using Base = Rp::Mask< D, RField<D>, FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::allocateBasis;
      using Base::allocateRGrid;
      using Base::setBasis;
      using Base::setRGrid;
      using Base::readBasis;
      using Base::readRGrid;
      using Base::basis;
      using Base::rgrid;
      using Base::phiTot;
      using Base::isAllocatedBasis;
      using Base::isAllocatedRGrid;
      using Base::hasData;
      using Base::isSymmetric;

   protected:

      /**
      * Calculate the average value of the rgrid_ member.
      */
      double rGridAverage() const override;

   };

   // Explicit instantiation declaration
   extern template class Mask<1>;
   extern template class Mask<2>;
   extern template class Mask<3>;

} // namespace Rpc

namespace Rp {
   // Explicit instantiation declaration for base class
   extern template 
   class Mask< 1, Prdc::Cpu::RField<1>, Rpc::FieldIo<1> >;
   extern template 
   class Mask< 2, Prdc::Cpu::RField<2>, Rpc::FieldIo<2> >;
   extern template 
   class Mask< 3, Prdc::Cpu::RField<3>, Rpc::FieldIo<3> >;
}

} // namespace Pscf
#endif

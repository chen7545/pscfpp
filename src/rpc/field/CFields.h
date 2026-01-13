#ifndef RPC_C_FIELDS_H
#define RPC_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>       // base class template
#include <rpc/field/FieldIo.h>      // base class template argument
#include <prdc/cpu/RField.h>        // base class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A list of c fields stored in both basis and r-grid format.
   *
   * This class is simply a named instantiation of the base class
   * template Rp::CFields, designed for use on CPU hardware. Please
   * see documentation of this base class for API documentation.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class CFields : public Rp::CFields<D, RField<D>, FieldIo<D> >
   {

   public:

      /// Alias for direct base class.
      using Base = Rp::CFields<D, RField<D>, FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::setNMonomer;
      using Base::setWriteUnitCell;
      using Base::allocateRGrid;
      using Base::allocateBasis;
      using Base::allocate;
      using Base::basis;
      using Base::rgrid;
      using Base::isAllocatedRGrid;
      using Base::isAllocatedBasis;
      using Base::hasData;
      using Base::isSymmetric;
      using Base::setHasData;
      using Base::setIsSymmetric;

   };

   // Explicit instantiation declarations
   extern template class CFields<1>;
   extern template class CFields<2>;
   extern template class CFields<3>;

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations of base class
namespace Pscf {
   namespace Rp {
      extern template 
      class CFields<1, Prdc::Cpu::RField<1>, Rpc::FieldIo<1> >;
      extern template 
      class CFields<2, Prdc::Cpu::RField<2>, Rpc::FieldIo<2> >;
      extern template 
      class CFields<3, Prdc::Cpu::RField<3>, Rpc::FieldIo<3> >;
   }
} 
#endif

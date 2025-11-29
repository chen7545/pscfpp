#ifndef CPC_C_FIELDS_H
#define CPC_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cl/CFields.h>     // base class template
#include <prdc/cpu/CField.h>     // template parameter
#include <cpc/field/FieldIo.h>   // template parameter

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A container of w fields.
   *
   * The public interface of this class is identical to that of the base
   * class template Pscf::Cp::CFields. Please see documentation of
   * that base class for API documentation.
   *
   * \ingroup Cpc_Field_Module
   */
   template <int D>
   class CFields 
     : public Cp::CFields<D, Prdc::Cpu::CField<D>, Cpc::FieldIo<D> >
   {
   public:

      /// Alias for base class.
      using Base = Cp::CFields< D, Prdc::Cpu::CField<D>, Cpc::FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::allocate;
      using Base::fields;
      using Base::field;
      using Base::writeFields;
      using Base::isAllocated;
      using Base::hasData;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nMonomer;
      using Base::fieldIo;

   };

   // Explicit instantiation declarations
   extern template class CFields<1>;
   extern template class CFields<2>;
   extern template class CFields<3>;

} // namespace Cpc
} // namespace Pscf

namespace Pscf {
namespace Cp {
   // Explicit instantiation declarations for base class
   extern template 
   class CFields<1, Prdc::Cpu::CField<1>, Cpc::FieldIo<1> >;
   extern template 
   class CFields<2, Prdc::Cpu::CField<2>, Cpc::FieldIo<2> >;
   extern template 
   class CFields<3, Prdc::Cpu::CField<3>, Cpc::FieldIo<3> >;
} 
} // namespace Pscf
#endif

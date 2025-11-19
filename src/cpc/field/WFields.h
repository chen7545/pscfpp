#ifndef CPC_W_FIELD_CONTAINER_H
#define CPC_W_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cl/WFields.h>     // base class template
#include <prdc/cpu/CField.h>     // template parameter
#include <cpc/field/FieldIo.h>   // template parameter

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * A container of w fields.
   *
   * The public interface of this class is identical to that of the base
   * class template Pscf::Prdc::Cl::WFields. Please see documentation
   * of that base class for API documentation.
   *
   * \ingroup Cpc_Field_Module
   */
   template <int D>
   class WFields 
     : public Cl::WFields<D, Prdc::Cpu::CField<D>, Cpc::FieldIo<D> >
   {
   public:

      /// Alias for base class.
      using Base = Cl::WFields< D, Prdc::Cpu::CField<D>, Cpc::FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::allocate;
      using Base::setFields;
      using Base::clear;
      using Base::operator();
      using Base::isAllocated;
      using Base::hasData;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nMonomer;
      using Base::fieldIo;

   private:

      /**
      * Assign one CField<D> to another: lhs = rhs.
      *
      * \param lhs  left-hand side of assignment
      * \param rhs  right-hand side of assignment
      */
      void assignField(CField<D>& lhs, CField<D> const & rhs) const 
      override;

   };

   // Explicit instantiation declarations
   extern template class WFields<1>;
   extern template class WFields<2>;
   extern template class WFields<3>;

} // namespace Cpc
} // namespace Pscf

namespace Pscf {
namespace Prdc {
   // Explicit instantiation declarations for base class
   extern template 
   class Cl::WFields<1, Cpu::CField<1>, Cpc::FieldIo<1> >;
   extern template 
   class Cl::WFields<2, Cpu::CField<2>, Cpc::FieldIo<2> >;
   extern template 
   class Cl::WFields<3, Cpu::CField<3>, Cpc::FieldIo<3> >;
} 
} // namespace Pscf
#endif

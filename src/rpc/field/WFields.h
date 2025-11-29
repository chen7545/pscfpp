#ifndef RPC_W_FIELD_CONTAINER_H
#define RPC_W_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/WFields.h>          // base class template
#include <prdc/cpu/RField.h>     // base class template argument
#include <rpc/field/FieldIo.h>   // base class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * The public interface of this class is identical to that of the base
   * class template Rp::WFields. Please see documentation of that base 
   * class template for API documentation.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class WFields 
     : public Rp::WFields<D, RField<D>, FieldIo<D> >
   {
   public:

      /// Alias for base class.
      using Base = Rp::WFields< D, RField<D>, FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::setNMonomer;
      using Base::allocateRGrid;
      using Base::allocateBasis;
      using Base::allocate;
      using Base::setBasis;
      using Base::setRGrid;
      using Base::readBasis;
      using Base::readRGrid;
      using Base::symmetrize;
      using Base::clear;
      using Base::basis;
      using Base::rgrid;
      using Base::isAllocatedRGrid;
      using Base::isAllocatedBasis;
      using Base::hasData;
      using Base::isSymmetric;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nBasis;
      using Base::nMonomer;
      using Base::fieldIo;

   private:

      /**
      * Assign one RField<D> to another: lhs = rhs.
      *
      * \param lhs  left-hand side of assignment
      * \param rhs  right-hand side of assignment
      */
      void assignRField(RField<D>& lhs, RField<D> const & rhs) const 
      override;

   };

   // Explicit instantiation declarations
   extern template class WFields<1>;
   extern template class WFields<2>;
   extern template class WFields<3>;

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations for base class
namespace Pscf {
   namespace Rp {
      extern template 
      class WFields<1, Prdc::Cpu::RField<1>, Rpc::FieldIo<1> >;
      extern template 
      class WFields<2, Prdc::Cpu::RField<2>, Rpc::FieldIo<2> >;
      extern template 
      class WFields<3, Prdc::Cpu::RField<3>, Rpc::FieldIo<3> >;
   }
} // namespace Pscf
#endif

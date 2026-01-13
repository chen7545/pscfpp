#ifndef RPG_MASK_H
#define RPG_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Mask.h>       // base class template
#include "FieldIo.h"             // base class template argument
#include <prdc/cuda/RField.h>    // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * A field to which the total monomer concentration is constrained.
   *
   * Please refer to the documentation of the Rp::Mask base class
   * template for more complete API documentation for this class template.
   * The public interface of Rpg::Mask is identical to that of the base
   * class template Rp::Mask.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Mask
     : public Rp::Mask< D, RField<D>, FieldIo<D> >
   {

   public:

      /// Base class typedef
      using Base = Rp::Mask< D, Prdc::Cuda::RField<D>, FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::allocateBasis;
      using Base::allocateRGrid;
      using Base::setBasis;
      using Base::setRGrid;
      using Base::readBasis;
      using Base::readRGrid;
      using Base::writeBasis;
      using Base::writeRGrid;
      using Base::basis;
      using Base::rgrid;
      using Base::phiTot;
      using Base::isAllocatedBasis;
      using Base::isAllocatedRGrid;
      using Base::hasData;
      using Base::isSymmetric;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nBasis;
      using Base::fieldIo;

      /**
      * Calculate the average value of the rgrid_ member.
      */
      double rGridAverage() const override;

   };

   // Explicit instantiation declarations
   extern template class Mask<1>;
   extern template class Mask<2>;
   extern template class Mask<3>;

} // namespace Rpg

namespace Rp {
   // Explicit instantiation declarations for base class template
   using namespace Prdc::Cuda;
   extern template class Mask< 1, RField<1>, Rpg::FieldIo<1> >;
   extern template class Mask< 2, RField<2>, Rpg::FieldIo<2> >;
   extern template class Mask< 3, RField<3>, Rpg::FieldIo<3> >;
}

} // namespace Pscf
#endif

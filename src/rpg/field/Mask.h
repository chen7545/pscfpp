#ifndef RPG_MASK_H
#define RPG_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/rl/Mask.h>        // base class template
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
   * Please refer to the documentation of the Rl::Mask base class
   * template for more complete API documentation for this class template.
   * The public interface of Rpg::Mask is identical to that of the base
   * class template Rl::Mask.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Mask
     : public Rl::Mask< D, RField<D>, FieldIo<D> >
   {

   public:

      /// Base class typedef
      using Base = Rl::Mask< D, Prdc::Cuda::RField<D>, FieldIo<D> >;

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

namespace Prdc {
   // Explicit instantiation declarations for base class template
   extern template class Rl::Mask< 1, Cuda::RField<1>, Rpg::FieldIo<1> >;
   extern template class Rl::Mask< 2, Cuda::RField<2>, Rpg::FieldIo<2> >;
   extern template class Rl::Mask< 3, Cuda::RField<3>, Rpg::FieldIo<3> >;
}

} // namespace Pscf
#endif

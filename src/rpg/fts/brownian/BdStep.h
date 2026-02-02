#ifndef RPG_BD_STEP_H
#define RPG_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>
#include <rpg/system/Types.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class BdStep : public Rp::BdStep<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator<D> object
      */
      BdStep(BdSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdStep< 1, Rpg::Types<1> >;
      extern template class BdStep< 2, Rpg::Types<2> >;
      extern template class BdStep< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BdStep<1>;
      extern template class BdStep<2>;
      extern template class BdStep<3>;
   }
}
#endif

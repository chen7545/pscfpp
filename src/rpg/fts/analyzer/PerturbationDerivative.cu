/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rp/fts/analyzer/PerturbationDerivative.tpp>

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   PerturbationDerivative<D>::PerturbationDerivative(
                                     Simulator<D>& simulator,
                                     System<D>& system)
    : Rp::PerturbationDerivative< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationDerivative< 1, Rpg::Types<1> >;
      template class PerturbationDerivative< 2, Rpg::Types<2> >;
      template class PerturbationDerivative< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class PerturbationDerivative<1>;
      template class PerturbationDerivative<2>;
      template class PerturbationDerivative<3>;
   }
}

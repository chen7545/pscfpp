/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimState.h"
#include <rp/fts/simulator/SimState.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template struct SimState<1, Prdc::Cuda::RField<1> >;
      template struct SimState<2, Prdc::Cuda::RField<2> >;
      template struct SimState<3, Prdc::Cuda::RField<3> >;
   }
   namespace Rpg {
      template struct SimState<1>;
      template struct SimState<2>;
      template struct SimState<3>;
   }
}

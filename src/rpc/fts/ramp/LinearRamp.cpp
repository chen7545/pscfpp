/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearRamp.h"
#include <rpc/fts/simulator/Simulator.h>
#include <util/global.h>

#if 0
#include <rpc/system/System.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/solvers/Block.h>
#include <rpc/field/Domain.h>
#endif

#include <rp/fts/ramp/LinearRamp.tpp>    // base class implementation

namespace Pscf {
namespace Rpc {

   // Constructor - creates association with a Simulator.
   template <int D>
   LinearRamp<D>::LinearRamp(Simulator<D>& simulator)
    : Rp::LinearRamp<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LinearRamp<1, Rpc::Types<1> >;
      template class LinearRamp<2, Rpc::Types<2> >;
      template class LinearRamp<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class LinearRamp<1>;
      template class LinearRamp<2>;
      template class LinearRamp<3>;
   }
}

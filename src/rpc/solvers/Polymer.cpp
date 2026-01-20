/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rpc/solvers/Block.h>
#include <rpc/solvers/Propagator.h>
#include <prdc/cpu/RField.h>
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>

#include <rp/solvers/Polymer.tpp>

namespace Pscf {
   namespace Rp {
      template class Polymer< 1, Rpc::Types<1> >;
      template class Polymer< 2, Rpc::Types<2> >;
      template class Polymer< 3, Rpc::Types<3> >;
   }

   namespace Rpc {
      template class Polymer<1>;
      template class Polymer<2>;
      template class Polymer<3>;
   }

}

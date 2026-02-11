/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"                    // header
#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rp/fts/analyzer/HamiltonianAnalyzer.tpp>  // implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(Simulator<D>& simulator,
                                   System<D>& system)
    : Rp::HamiltonianAnalyzer< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class HamiltonianAnalyzer< 1, Rpc::Types<1> >;
      template class HamiltonianAnalyzer< 2, Rpc::Types<2> >;
      template class HamiltonianAnalyzer< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class HamiltonianAnalyzer<1>;
      template class HamiltonianAnalyzer<2>;
      template class HamiltonianAnalyzer<3>;
   }
}

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.tpp"

namespace Pscf {
   namespace Rp {
      template class BdSimulator<1, Rpc::Types<1> >;
      template class BdSimulator<2, Rpc::Types<2> >;
      template class BdSimulator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class BdSimulator<1>;
      template class BdSimulator<2>;
      template class BdSimulator<3>;
   }
}

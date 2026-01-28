/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.tpp"

namespace Pscf { 
   namespace Rp {
      template class Propagator<1, Rpc::Types<1> >;
      template class Propagator<2, Rpc::Types<2> >;
      template class Propagator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class Propagator<1>;
      template class Propagator<2>;
      template class Propagator<3>;
   }
}

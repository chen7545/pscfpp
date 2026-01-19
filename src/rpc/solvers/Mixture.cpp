/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"

namespace Pscf {
   namespace Rp { 
      template class Mixture<1, Rpc::Types<1> >;
      template class Mixture<2, Rpc::Types<2> >;
      template class Mixture<3, Rpc::Types<3> >;
   }
   namespace Rpc { 
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}

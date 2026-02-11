/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.tpp"

namespace Pscf {
   namespace Rp {
      template class Simulator<1, Rpg::Types<1> >;
      template class Simulator<2, Rpg::Types<2> >;
      template class Simulator<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Simulator<1>;
      template class Simulator<2>;
      template class Simulator<3>;
   }
}

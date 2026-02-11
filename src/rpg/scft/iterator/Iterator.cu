/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.tpp"

namespace Pscf {
   namespace Rp {
      template class Iterator<1, Rpg::System<1> >;
      template class Iterator<2, Rpg::System<2> >;
      template class Iterator<3, Rpg::System<3> >;
   }
   namespace Rpg {
      template class Iterator<1>;
      template class Iterator<2>;
      template class Iterator<3>;
   }
}

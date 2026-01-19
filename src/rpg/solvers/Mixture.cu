/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"

namespace Pscf {
   namespace Rp { 
      // Explicit instantiation definitions for base class
      template class Mixture<1, Rpg::Types<1> >;
      template class Mixture<2, Rpg::Types<2> >;
      template class Mixture<3, Rpg::Types<3> >;
   }
   namespace Rpg { 
      // Explicit instantiation definitions for this class
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}

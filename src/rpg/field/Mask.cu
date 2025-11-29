/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.tpp"
#include <prdc/rl/Mask.tpp> // needed for implementation

namespace Pscf {
   namespace Rp {
      // Explicit instantiation of base class
      using namespace Prdc::Cuda;
      template class Mask< 1, RField<1>, Rpg::FieldIo<1> >;
      template class Mask< 2, RField<2>, Rpg::FieldIo<2> >;
      template class Mask< 3, RField<3>, Rpg::FieldIo<3> >;
   } 
   namespace Rpg {
      // Explicit instantiation of this class
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
      
   } 
} 

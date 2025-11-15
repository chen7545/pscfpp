/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"
#include <rpg/field/FieldIo.h>
#include <prdc/cuda/RField.h>
#include <prdc/rl/CFields.tpp>  // Base class implementation

namespace Pscf {
   namespace Prdc {
      // Explicit instantiations of base class template
      template class Rl::CFields<1, Cuda::RField<1>, Rpg::FieldIo<1> >;
      template class Rl::CFields<2, Cuda::RField<2>, Rpg::FieldIo<2> >;
      template class Rl::CFields<3, Cuda::RField<3>, Rpg::FieldIo<3> >;
   } 
   namespace Rpg {
      // Explicit instantiations of this class
      template class CFields<1>;
      template class CFields<2>;
      template class CFields<3>;
   } 
}

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"
#include <prdc/rl/CFields.tpp>
#include <prdc/cpu/RField.h>
#include <rpc/field/FieldIo.h>

namespace Pscf {
   namespace Rl {
      // Explicit instantiations of base class
      template class CFields<1, Prdc::Cpu::RField<1>, Rpc::FieldIo<1> >;
      template class CFields<2, Prdc::Cpu::RField<2>, Rpc::FieldIo<2> >;
      template class CFields<3, Prdc::Cpu::RField<3>, Rpc::FieldIo<3> >;
   }
   namespace Rpc {
      // Explicit instantiations of this class
      template class CFields<1>;
      template class CFields<2>;
      template class CFields<3>;

   }
}

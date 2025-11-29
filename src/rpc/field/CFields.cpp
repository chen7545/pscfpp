/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"              // class header
#include <rp/CFields.tpp>         // base class template implementation
#include <prdc/cpu/RField.h>      // base class template argument
#include <rpc/field/FieldIo.h>    // base class template argument

namespace Pscf {
   namespace Rp {
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

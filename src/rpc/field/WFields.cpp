/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.tpp"

namespace Pscf {
   namespace Rp {
      // Explicit instantiation definitions for base class
      using namespace Prdc::Cpu;
      template class Rp::WFields<1, RField<1>, Rpc::FieldIo<1> >;
      template class Rp::WFields<2, RField<2>, Rpc::FieldIo<2> >;
      template class Rp::WFields<3, RField<3>, Rpc::FieldIo<3> >;
   }
   namespace Rpc {
      // Explicit instantiation definitions
      template class WFields<1>;
      template class WFields<2>;
      template class WFields<3>;
   }
}

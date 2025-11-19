/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.tpp"

namespace Pscf {
namespace Prdc {
   // Explicit instantiation definitions for base class
   template class Rl::WFields<1, RField<1>, Rpc::FieldIo<1> >;
   template class Rl::WFields<2, RField<2>, Rpc::FieldIo<2> >;
   template class Rl::WFields<3, RField<3>, Rpc::FieldIo<3> >;
}
namespace Rpc {
   // Explicit instantiation definitions
   template class WFields<1>;
   template class WFields<2>;
   template class WFields<3>;
}
}

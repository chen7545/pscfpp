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
   template class Cl::WFields<1, Cpu::CField<1>, Cpc::FieldIo<1> >;
   template class Cl::WFields<2, Cpu::CField<2>, Cpc::FieldIo<2> >;
   template class Cl::WFields<3, Cpu::CField<3>, Cpc::FieldIo<3> >;
}
namespace Cpc {
   // Explicit instantiation definitions
   template class WFields<1>;
   template class WFields<2>;
   template class WFields<3>;
}
}

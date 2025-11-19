/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"
#include <prdc/cl/CFields.tpp>
#include <prdc/cpu/CField.h>
#include <cpc/field/FieldIo.h>

namespace Pscf {
namespace Prdc {
   // Explicit instantiation definitions for base class
   template class Cl::CFields<1, Cpu::CField<1>, Cpc::FieldIo<1> >;
   template class Cl::CFields<2, Cpu::CField<2>, Cpc::FieldIo<2> >;
   template class Cl::CFields<3, Cpu::CField<3>, Cpc::FieldIo<3> >;
}
namespace Cpc {
   // Explicit instantiation definitions
   template class CFields<1>;
   template class CFields<2>;
   template class CFields<3>;
}
}

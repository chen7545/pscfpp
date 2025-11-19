/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"

namespace Pscf {
   namespace Prdc {
      template class Cl::FieldIo<1, Cpu::CField<1>, Cpu::FFT<1> >;
      template class Cl::FieldIo<2, Cpu::CField<2>, Cpu::FFT<2> >;
      template class Cl::FieldIo<3, Cpu::CField<3>, Cpu::FFT<3> >;
   }
   namespace Cpc {
      template class FieldIo<1>;
      template class FieldIo<2>;
      template class FieldIo<3>;
   } 
}

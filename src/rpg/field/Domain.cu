/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.tpp"

namespace Pscf {
   namespace Prdc {
      using namespace Cuda;
      template class Rl::Domain<1, FFT<1>, WaveList<1>, Rpg::FieldIo<1> >;
      template class Rl::Domain<2, FFT<2>, WaveList<2>, Rpg::FieldIo<2> >;
      template class Rl::Domain<3, FFT<3>, WaveList<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      template class Domain<1>;
      template class Domain<2>;
      template class Domain<3>;
   }
}

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"
#include <util/math/Constants.h>

namespace Pscf {

   namespace Rp {
      // Explicit instantiation
      using namespace Prdc::Cuda;
      template class FieldIo<1, RField<1>, RFieldDft<1>, FFT<1> >;
      template class FieldIo<2, RField<2>, RFieldDft<2>, FFT<2> >;
      template class FieldIo<3, RField<3>, RFieldDft<3>, FFT<3> >;
   }
   
   namespace Rpg {
      // Explicit instantiations
      template class FieldIo<1>;
      template class FieldIo<2>;
      template class FieldIo<3>;
   } 

} // namespace Pscf

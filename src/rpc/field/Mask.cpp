/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.tpp"                 // class implementation
#include <prdc/rl/Mask.tpp>  // base class implementation

namespace Pscf {
   namespace Rp {
      // Explicit instantiation definitions for base class
      using namespace Prdc::Cpu;
      template class Rp::Mask< 1, RField<1>, Rpc::FieldIo<1> >;
      template class Rp::Mask< 2, RField<2>, Rpc::FieldIo<2> >;
      template class Rp::Mask< 3, RField<3>, Rpc::FieldIo<3> >;
   
   } // namespace Prdc
   namespace Rpc {
      // Explicit instantiation definitions for base class
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
   } // namespace Rpc
} // namespace Pscf

#ifndef RPC_W_FIELDS_TPP
#define RPC_W_FIELDS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"               // class header
#include <rpc/field/FieldIo.tpp>
#include <rp/field/WFields.tpp>    // base class implementation

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   template <int D>
   void WFields<D>::assignRField(RField<D>& lhs,
                                 RField<D> const & rhs) const
   {
      int n = rhs.capacity();
      UTIL_CHECK(lhs.capacity() == n);
      UTIL_CHECK(RpWFields::meshSize() == n);
      for (int i = 0; i < n; ++i) {
         lhs[i] = rhs[i];
      }
   }

} // namespace Rpc
} // namespace Pscf
#endif

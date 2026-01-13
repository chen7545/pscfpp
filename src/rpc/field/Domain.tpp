#ifndef RPC_DOMAIN_TPP
#define RPC_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"                 // class header
#include <rp/field/Domain.tpp>      // base class template implementation
#include <rpc/field/FieldIo.tpp>    // base class template argument
#include <prdc/cpu/WaveList.h>      // base class template argument
#include <prdc/cpu/FFT.h>           // base class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : Rp::Domain<D, FFT<D>, WaveList<D>, FieldIo<D> >()
   {  ParamComposite::setClassName("Domain"); }

} // namespace Rpc
} // namespace Pscf
#endif

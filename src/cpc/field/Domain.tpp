#ifndef CPC_DOMAIN_TPP
#define CPC_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <prdc/cl/Domain.tpp>

#include <cpc/field/FieldIo.tpp>
#include <prdc/cpu/WaveList.h>
#include <prdc/cpu/FFT.h>

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : Cl::Domain<D, FFT<D>, WaveList<D>, FieldIo<D> >()
   {  ParamComposite::setClassName("Domain"); }

} // namespace Cpc
} // namespace Pscf
#endif

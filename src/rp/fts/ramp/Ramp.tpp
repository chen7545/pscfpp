#ifndef RP_RAMP_TPP
#define RP_RAMP_TPP

#include "Ramp.h"
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   Ramp<D>::Ramp(typename T::Simulator& simulator)
    : ParamComposite(),
      simulatorPtr_(&simulator)
   {}
   
   /* 
   * Destructor.
   */
   template <int D>
   Ramp<D>::~Ramp()
   {}

   /* 
   * Setup before simulation - sets the nStep member variable. 
   */
   template <int D>
   void Ramp<D>::setup(int nStep)
   {  nStep_ = nStep; }

}
}
#endif 

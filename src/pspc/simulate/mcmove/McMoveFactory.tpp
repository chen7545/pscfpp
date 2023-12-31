#ifndef PSPC_MC_MOVE_FACTORY_TPP
#define PSPC_MC_MOVE_FACTORY_TPP

#include "McMoveFactory.h"  
#include <pspc/simulate/mcmove/McSimulator.h>

// Subclasses of McMove 
#include "RealMove.h"
#include "RealWcMove.h"
#include "FourierMove.h"
#include "FourierWcMove.h"
#include "ForceBiasMove.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   McMoveFactory<D>::McMoveFactory(McSimulator<D>& mcSimulator)
    : mcSimulatorPtr_(&mcSimulator)
   {}

   /* 
   * Return a pointer to a instance of McMove subclass className.
   */
   template <int D>
   McMove<D>* McMoveFactory<D>::factory(const std::string &className) const
   {
      McMove<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      
      // Try to match classname
      if (className == "RealMove") {
         ptr = new RealMove<D>(*mcSimulatorPtr_);
      } else if (className == "RealWcMove") {
         ptr = new RealWcMove<D>(*mcSimulatorPtr_);
      } else if (className == "FourierMove") {
         ptr = new FourierMove<D>(*mcSimulatorPtr_);
      } else if (className == "FourierWcMove") {
         ptr = new FourierWcMove<D>(*mcSimulatorPtr_);
      } else if (className == "ForceBiasMove") {
         ptr = new ForceBiasMove<D>(*mcSimulatorPtr_);
      }

      return ptr;
   }

}
}
#endif

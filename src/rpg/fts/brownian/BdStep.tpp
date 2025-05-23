#ifndef RPG_BD_STEP_TPP
#define RPG_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStep.h"

#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/System.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BdStep<D>::BdStep(BdSimulator<D>& simulator)
    : simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      cudaRandomPtr_(&(simulator.cudaRandom()))
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   BdStep<D>::~BdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void BdStep<D>::readParameters(std::istream &in)
   {}

   /*
   * Setup at beginning of loop.
   */
   template <int D>
   void BdStep<D>::setup()
   {}

   template <int D>
   void BdStep<D>::output()
   {}

   template<int D>
   void BdStep<D>::outputTimers(std::ostream& out)
   {}

   template<int D>
   void BdStep<D>::clearTimers()
   {}

}
}
#endif

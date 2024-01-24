#ifndef PSPC_MIDSTEP_BD_STEP_H
#define PSPC_MIDSTEP_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"

#include <prdc/cpu/RField.h>
#include <util/containers/DArray.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * Midstep predictor Brownian dynamics step.
   *
   * \ingroup Pspc_Simulate_BdStep_Module
   */
   template <int D>
   class MidstepBdStep : public BdStep<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param bdSimulator  parent BdSimulator object
      */
      MidstepBdStep(BdSimulator<D>& bdSimulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~MidstepBdStep();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Setup before simulation.
      */
      virtual void setup();

      /**
      * Take a single Brownian dynamics step.
      */
      virtual void step();

   protected:

      using BdStep<D>::system;
      using BdStep<D>::simulator;
      using BdStep<D>::random;
      using ParamComposite::read;

   private:

      // Local copies of w fields
      DArray< RField<D> > wh_;
      DArray< RField<D> > wf_;

      // Random displacement components
      DArray< RField<D> > eta_;

      // Change in one component of wc
      RField<D> dwc_;

      // Change in pressure field component 
      RField<D> dwp_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

   #ifndef PSPC_MIDSTEP_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class MidstepBdStep<1>;
   extern template class MidstepBdStep<2>;
   extern template class MidstepBdStep<3>;
   #endif

}
}
#endif

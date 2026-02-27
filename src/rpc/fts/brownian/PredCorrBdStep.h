#ifndef RPC_PRED_CORR_BD_STEP_H
#define RPC_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"
#include <prdc/cpu/RField.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Predictor-corrector Brownian dynamics stepper.
   *
   * \see \ref rp_PredCorrBdStep_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class PredCorrBdStep : public BdStep<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      PredCorrBdStep(BdSimulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~PredCorrBdStep();

      /**
      * Read body of parameter file block and initialize.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before simulation.
      */
      void setup() override;

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true if converged, false if failed to converge
      */
      bool step() override;

   protected:

      using BdStep<D>::system;
      using BdStep<D>::simulator;

   private:

      // Predicted values of fields (monomer fields)
      DArray< RField<D> > wp_;

      // Corrected (final) values of fields (monomer fields)
      DArray< RField<D> > wf_;

      // Initial deterministic forces (eigenvector components)
      DArray< RField<D> > dci_;

      // Random displacement components (eigenvector components)
      DArray< RField<D> > eta_;

      // Change in one component of wc
      RField<D> dwc_;

      // Change in pressure field component
      RField<D> dwp_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

   // Explicit instantiation declarations
   extern template class PredCorrBdStep<1>;
   extern template class PredCorrBdStep<2>;
   extern template class PredCorrBdStep<3>;

}
}
#endif

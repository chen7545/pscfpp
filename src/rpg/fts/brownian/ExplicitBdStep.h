#ifndef RPG_EXPLICIT_BD_STEP_H
#define RPG_EXPLICIT_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"                    // base class
#include <prdc/cuda/RField.h>          // member
#include <util/containers/DArray.h>    // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Explicit Euler-Mayurama Brownian dynamics step.
   *
   * \see \ref rpc_ExplicitBdStep_page "Manual Page"
   *
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class ExplicitBdStep : public BdStep<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      ExplicitBdStep(BdSimulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~ExplicitBdStep();

      /**
      * Read body of parameter file block. 
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before simulation.
      */
      void setup() override;

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true iff the compressor converged, false otherwise
      */
      bool step() override;

   protected:

      using BdStepT = BdStep<D>;
      using BdStep<D>::system;
      using BdStep<D>::simulator;
      using BdStep<D>::vecRandom;

   private:

      // Local copy of w fields
      DArray< RField<D> > w_;

      // Change in one component of wc
      RField<D> dwc_;

      /// Normal-distributed random fields
      RField<D> gaussianField_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

   // Explicit instantiation declarations
   extern template class ExplicitBdStep<1>;
   extern template class ExplicitBdStep<2>;
   extern template class ExplicitBdStep<3>;

}
}
#endif

#ifndef RPC_EXPLICIT_BD_STEP_H
#define RPC_EXPLICIT_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"                    // base class
#include <prdc/cpu/RField.h>           // member
#include <util/containers/DArray.h>    // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Explicit Euler-Mayurama Brownian dynamics step.
   *
   * \see \ref rpc_ExplicitBdStep_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Brownian_Module
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

      using BdStep<D>::system;
      using BdStep<D>::simulator;
      using BdStep<D>::random;
      using BdStep<D>::vecRandom;

   private:

      // Local copy of w fields
      DArray< RField<D> > w_;

      // Change in one component of wc
      RField<D> dwc_;

      // Normal distributed random numbers
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

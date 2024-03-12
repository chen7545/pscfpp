#ifndef RPC_EXPLICIT_BD_STEP_H
#define RPC_EXPLICIT_BD_STEP_H

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
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Rpc_Simulate_BdStep_Module
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
      *
      * Empty default implementation.
      */
      virtual ~ExplicitBdStep();

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
      
      /**
      * Return number of moves that fail to converge.
      */
      virtual long nFail();

   protected:

      using BdStep<D>::system;
      using BdStep<D>::simulator;
      using BdStep<D>::random;
      using ParamComposite::read;

   private:

      // Local copy of w fields
      DArray< RField<D> > w_;

      // Change in one component of wc
      RField<D> dwc_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;
      
      // Number of moves that fail to converge.
      long  nFail_;

   };
   
   /*
   * Return number of moves that fail to converge.
   */
   template <int D>
   inline long ExplicitBdStep<D>::nFail()
   {  return nFail_; }

   #ifndef RPC_EXPLICIT_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class ExplicitBdStep<1>;
   extern template class ExplicitBdStep<2>;
   extern template class ExplicitBdStep<3>;
   #endif

}
}
#endif

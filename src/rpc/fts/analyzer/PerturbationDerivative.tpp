#ifndef RPC_PERTURBATION_DERIVATIVE_TPP
#define RPC_PERTURBATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/perturbation/Perturbation.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   PerturbationDerivative<D>::PerturbationDerivative(Simulator<D>& simulator,
                                                     System<D>& system)
    : AverageAnalyzer<D>(simulator, system)
   {  ParamComposite::setClassName("PerturbationDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   PerturbationDerivative<D>::~PerturbationDerivative()
   {}

   template <int D>
   double PerturbationDerivative<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(simulator().hasPerturbation());

      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      return simulator().perturbation().df();
   }

   template <int D>
   void PerturbationDerivative<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         std::ofstream& outputFile = AverageAnalyzer<D>::outputFile_;
         UTIL_CHECK(outputFile.is_open());
         double lambda = simulator().perturbation().lambda();
         outputFile << Int(step);
         outputFile << Dbl(lambda);
         outputFile << Dbl(value);
         outputFile << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }

} // namespace Rpg
} // namespace Pscf
#endif

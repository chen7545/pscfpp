#ifndef RPC_CHI_DERIVATIVE_TPP
#define RPC_CHI_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/interaction/Interaction.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ChiDerivative<D>::ChiDerivative(Simulator<D>& simulator,
                                   System<D>& system)
    : AverageAnalyzer<D>(simulator, system)
   {  ParamComposite::setClassName("ChiDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   ChiDerivative<D>::~ChiDerivative()
   {}

   template <int D>
   double ChiDerivative<D>::compute()
   {
      // Preconditions
      UTIL_CHECK(system().w().hasData());
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);
      // Only valid for nMonomer == 2

      const double vSystem = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();
      double chi = system().interaction().chi(0,1);

      // Pre-requisite computations
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      // Get field Hamiltonian per monomer
      double hField = simulator().fieldHamiltonian();
      hField /= (double)nMonomerSystem;

      // Compute derivative of f (per monomer) w/respect to bare chi
      double dfdchi = -(hField - 0.5*simulator().sc(nMonomer - 1))/chi
                    + 0.25;

      // Convert to derivative of total free energy F
      dfdchi *= (double)nMonomerSystem;

      // With N term
      dfdchi += double(meshSize)/(2.0 * chi);

      return dfdchi;
   }

   template <int D>
   void ChiDerivative<D>::outputValue(int step, double value)
   {
      int nSamplePerOutput = AverageAnalyzer<D>::nSamplePerOutput();
      if (simulator().hasRamp() && nSamplePerOutput == 1) {
         UTIL_CHECK(outputFile_.is_open());
         double chi = system().interaction().chi(0,1);
         outputFile_ << Int(step);
         outputFile_ << Dbl(chi);
         outputFile_ << Dbl(value) << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }

}
}
#endif

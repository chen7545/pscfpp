#ifndef RPG_HAMILTONIAN_ANALYZER_TPP
#define RPG_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <iostream>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(Simulator<D>& simulator,
                                               System<D>& system)
    : AverageListAnalyzer<D>(simulator, system),
      idealId_(-1),
      fieldId_(-1),
      totalId_(-1)
   {  ParamComposite::setClassName("HamiltonianAnalyzer"); }

   /*
   * Read interval and outputFileName.
   */
   template <int D>
   void HamiltonianAnalyzer<D>::readParameters(std::istream& in)
   {
      AverageListAnalyzer<D>::readParameters(in);
      AverageListAnalyzer<D>::initializeAccumulators(3);

      idealId_ = 0;
      AverageListAnalyzer<D>::setName(idealId_, "ideal");
      fieldId_ = 1;
      AverageListAnalyzer<D>::setName(fieldId_, "field");
      totalId_ = 2;
      AverageListAnalyzer<D>::setName(totalId_, "total");
   }

   /*
   * Output energy to file.
   */
   template <int D>
   void HamiltonianAnalyzer<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());
      if (!system().c().hasData()) {
         system().compute();
      }
      UTIL_CHECK(system().c().hasData());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      UTIL_CHECK(simulator().hasWc());
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      UTIL_CHECK(simulator().hasHamiltonian());

      double ideal = simulator().idealHamiltonian();
      AverageListAnalyzer<D>::setValue(idealId_, ideal);

      double field = simulator().fieldHamiltonian();
      AverageListAnalyzer<D>::setValue(fieldId_, field);

      double total = simulator().hamiltonian();
      AverageListAnalyzer<D>::setValue(totalId_, total);
   }

}
}
#endif

#ifndef RPG_BOX_LENGTH_DERIVATIVE_TPP
#define RPG_BOX_LENGTH_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BoxLengthDerivative.h"

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/cuda/resources.h>
#include <util/param/ParamComposite.h>  

namespace Pscf {
namespace Rpg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BoxLengthDerivative<D>::BoxLengthDerivative(Simulator<D>& simulator, 
                                   System<D>& system) 
    : AverageAnalyzer<D>(simulator, system)
   {  setClassName("BoxLengthDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   BoxLengthDerivative<D>::~BoxLengthDerivative() 
   {}

   template <int D>
   double BoxLengthDerivative<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());
      
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      // Simulations for cubic box
      UTIL_CHECK(D == 3);

      int nParameter = system().domain().unitCell().nParameter();
      UTIL_CHECK(nParameter == 1);
      
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();

      /* 
      * Compute field Hamiltonian per monomer.
      * The fieldHamitonian is calculated in the computeHamiltonian() function,
      * located in rpg/fts/Simulator.tpp 
      */
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      // Box length
      double l = system().domain().unitCell().parameter(0);

      // Obtain fieldHamiltonian 
      double HW = simulator().fieldHamiltonian();

      // The fieldHamiltonian contribution to the derivative
      double dFdL = HW/vSystem * 3.0 * l * l;

      // Obtain ideal gas Hamiltonian - n[lnQ_id + W_{+}/M]
      double nlnQ= simulator().idealHamiltonian();
      dFdL += 3.0/l * nlnQ;

      // Obtain stress -1/Q dQdl per monomer
      if (!system().mixture().hasStress()) {
         system().computeStress();
      }
      double stress = system().mixture().stress(0);
      dFdL += stress * nMonomerSystem;

      // With N term
      dFdL -= 3.0 * double(meshSize)/(2.0 * l);

      return dFdL;
   }
   
   template <int D>
   void BoxLengthDerivative<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double l = system().domain().unitCell().parameter(0);
         
         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(l);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }
   
}
}
#endif

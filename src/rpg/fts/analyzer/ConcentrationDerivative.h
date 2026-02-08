#ifndef RPG_CONCENTRATION_DERIVATIVE_H
#define RPG_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"         // base class template

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to concentration.
   *
   * \see \ref \rp_ConcentrationDerivative_page "Manual Page"
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationDerivative : public AverageAnalyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      ConcentrationDerivative(Simulator<D>& simulator, System<D>& system);

   protected:

      /**
      * Compute and return the derivative of H w/ respect to concentration.
      */
      virtual double compute();

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);

      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;

   };

   // Explicit instantiation declarations
   extern template class ConcentrationDerivative<1>;
   extern template class ConcentrationDerivative<2>;
   extern template class ConcentrationDerivative<3>;

}
}
#endif

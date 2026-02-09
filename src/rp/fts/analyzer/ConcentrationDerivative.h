#ifndef RP_CONCENTRATION_DERIVATIVE_H
#define RP_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to concentration.
   *
   * \see \ref rp_ConcentrationDerivative_page "Manual Page"
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class ConcentrationDerivative : public T::AverageAnalyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationDerivative(typename T::Simulator& simulator,
                              typename T::System& system);

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

      using AnayzerT = typename T::Analyzer;
      using AverageAnayzerT = typename T::AverageAnalyzer;
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   };

} // namespace Rp
} // namespace Pscf
#endif

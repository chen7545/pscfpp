#ifndef RP_PERTURBATION_DERIVATIVE_H
#define RP_PERTURBATION_DERIVATIVE_H

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
   * Evaluate derivative of H w/ respect to perturbation parameter lambda.
   *
   * \see rp_PerturbationDerivative_page "Manual Page"
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class PerturbationDerivative : public T::AverageAnalyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      PerturbationDerivative(typename T::Simulator& simulator, 
                             typename T::System& system);

      /**
      * Destructor.
      */
      virtual ~PerturbationDerivative();

   protected:

      /**
      * Compute and return the derivative of H w/ respect to lambda.
      */
      virtual double compute();

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);

      using AnalyzerT = typename T::Analyzer;
      using AverageAnalyzerT = typename T::AverageAnalyzer;
      using AverageT::simulator;
      using AverageT::system;

   };

} // namespace Rp
} // namespace Pscf
#endif

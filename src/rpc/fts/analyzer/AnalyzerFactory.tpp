#ifndef RPC_ANALYZER_FACTORY_TPP
#define RPC_ANALYZER_FACTORY_TPP

#include "AnalyzerFactory.h"

// Subclasses of Analyzer
#include "TrajectoryWriter.h"
#include "ConcentrationWriter.h"
#include "HamiltonianAnalyzer.h"
#include "BinaryStructureFactorGrid.h"
#include "StepLogger.h"
#include "PerturbationDerivative.h"
#include "ChiDerivative.h"
#include "ConcentrationDerivative.h"
#include "MaxOrderParameter.h"
#include "FourthOrderParameter.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AnalyzerFactory<D>::AnalyzerFactory(Simulator<D>& simulator,
                                       System<D>& system)
    : simPtr_(&simulator),
      sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Analyzer subclass className.
   */
   template <int D>
   Analyzer<D>* AnalyzerFactory<D>::factory(const std::string &className) const
   {
      Analyzer<D>* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;


      // Try to match classname
      if (className == "TrajectoryWriter") {
         ptr = new TrajectoryWriter<D>(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationWriter") {
         ptr = new ConcentrationWriter<D>(*simPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new HamiltonianAnalyzer<D>(*simPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactorGrid") {
         ptr = new BinaryStructureFactorGrid<D>(*simPtr_, *sysPtr_);
      } else if (className == "StepLogger") {
         ptr = new StepLogger<D>(*simPtr_, *sysPtr_);
      } else if (className == "PerturbationDerivative") {
         ptr = new PerturbationDerivative<D>(*simPtr_, *sysPtr_);
      } else if (className == "ChiDerivative") {
         ptr = new ChiDerivative<D>(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationDerivative") {
         ptr = new ConcentrationDerivative<D>(*simPtr_, *sysPtr_);
      } else if (className == "MaxOrderParameter") {
         ptr = new MaxOrderParameter<D>(*simPtr_, *sysPtr_);
      } else if (className == "FourthOrderParameter") {
         ptr = new FourthOrderParameter<D>(*simPtr_, *sysPtr_);
      }

      return ptr;
   }

}
}
#endif

#ifndef RP_ANALYZER_MANAGER_TPP
#define RP_ANALYZER_MANAGER_TPP

#include "AnalyzerManager.h"
#include <util/param/Factory.h>
#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   AnalyzerManager<D>::AnalyzerManager(typename T::Simulator& simulator,
                                       typename T::System& system)
    : Base(),
      simulatorPtr_(&simulator),
      systemPtr_(&system)
   {  ParamComposite::setClassName("AnalyzerManager"); }

   /*
   * Destructor.
   */
   template <int D>
   AnalyzerManager<D>::~AnalyzerManager()
   {}

   /*
   * Return a pointer to a new AnalyzerFactory object.
   */
   template <int D>
   Factory<typename T::Analyzer>* AnalyzerManager<D>::newDefaultFactory() 
   const
   {  return new AnalyzerFactory<D>(*simulatorPtr_, *systemPtr_); }

   /*
   * Read body of parameter file block.
   */
   template <int D>
   void AnalyzerManager<D>::readParameters(std::istream &in)
   {
      AnalyzerT::baseInterval = 1;
      ParamComposite::readOptional(in, "baseInterval",
                                   AnalyzerT::baseInterval);
      Base::readParameters(in);
   }

   /*
   * Call setup method of each analyzer.
   */
   template <int D>
   void AnalyzerManager<D>::setup()
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].setup();
      }
   }

   /*
   * Call sample method of each analyzer.
   */
   template <int D>
   void AnalyzerManager<D>::sample(long iStep)
   {
      int baseInterval = AnalyzerT::baseInterval;
      UTIL_CHECK(baseInterval > 0);
      UTIL_CHECK(iStep % baseInterval == 0);
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].sample(iStep);
      }
   }

   /*
   * Call output method of each analyzer.
   */
   template <int D>
   void AnalyzerManager<D>::output()
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].output();
      }
   }

}
}
#endif

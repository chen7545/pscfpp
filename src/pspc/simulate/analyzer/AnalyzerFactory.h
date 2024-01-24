#ifndef PSPC_ANALYZER_FACTORY_H
#define PSPC_ANALYZER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/simulate/analyzer/Analyzer.h>
#include <util/param/Factory.h>  
#include <string>

namespace Pscf {
namespace Pspc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Factory for subclasses of Analyzer.
   *
   * \ingroup Pspc_Simulate_Analyzer_Module
   */
   template <int D>
   class AnalyzerFactory : public Factory< Analyzer<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator<D> object
      * \param system  parent System<D> object
      */
      AnalyzerFactory(Simulator<D>& simulator, System<D>& system);

      /**
      * Method to create any Analyzer supplied with PSCF.
      *
      * \param className name of the Analyzer subclass
      * \return Analyzer* pointer to new instance of className
      */
      Analyzer<D>* factory(const std::string &className) const;

      using Factory< Analyzer<D> >::trySubfactories;

   private:
      
      /// Pointer to the parent system.
      System<D>* sysPtr_;
      
      /// Pointer to the parent simulator.
      Simulator<D>* simulatorPtr_;

   };

   #ifndef PSPC_ANALYZER_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class AnalyzerFactory<1>;
   extern template class AnalyzerFactory<2>;
   extern template class AnalyzerFactory<3>;
   #endif

}
}
#endif

#ifndef RPC_ANALYZER_MANAGER_H
#define RPC_ANALYZER_MANAGER_H

#include "Analyzer.h"                  // template parameter
#include <util/param/Manager.h>        // base class template

// Forward declarations
namespace Util {
   template <class T> class Factory;
}
namespace Pscf {
   namespace Rpc {
      template <int D> class System;
      template <int D> class Simulator;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Manager for a list of Analyzer objects.
   *
   * \see \ref rp_AnalyzerManager_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AnalyzerManager : public Manager< Analyzer<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator
      * \param system  parent System
      */
      AnalyzerManager(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~AnalyzerManager();

      /**
      * Read body of parameter file block.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Call the setup function of each Analyzer.
      *
      * This function should be called just before the main
      * simulation loop, after an initial configuration is
      * known. It calls the setup() functionfor each
      * analyzer, or does nothing if size() == 0.
      */
      void setup();

      /**
      * Call the sample function of each Analyzer.
      *
      * \pre Analyzer<D>::baseInterval > 0
      * \pre iStep % Analyzer<D>::baseInterval == 0
      *
      * \param iStep  step counter for main loop
      */
      void sample(long iStep);

      /**
      * Call the output function of each analyzer.
      *
      * This function should be called after the main
      * simulation loop. It calls the output() function
      * of each analyzer, or does nothing if size() == 0.
      */
      void output();

   private:

      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;

      /**
      * Pointer to parent System.
      */
      System<D>* systemPtr_;

      /**
      * Return pointer to a new AnalyzerFactory.
      */
      Factory< Analyzer<D> >* newDefaultFactory() const override;

      using Base = Manager< Analyzer<D> >;

   };

   // Explicit instantiation declarations
   extern template class AnalyzerManager<1>;
   extern template class AnalyzerManager<2>;
   extern template class AnalyzerManager<3>;

}
}
#endif

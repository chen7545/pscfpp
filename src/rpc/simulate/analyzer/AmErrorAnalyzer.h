#ifndef RPC_AM_ERROR_ANALYZER_H
#define RPC_AM_ERROR_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/accumulators/Average.h>           // member
#include <rpc/System.h>
#include <rpc/simulate/Simulator.h>
#include <iostream>
#include <vector>

namespace Pscf {
namespace Rpc 
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of several sampled real variables, and
   * optionally writes block averages to a data file during a simulation. 
   * It is intended for use as a base class for Analyzers that evaluate 
   * averages and (optionally) block averages for specific physical variables.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class AmErrorAnalyzer : public Analyzer<D>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent SystemType object. 
      */
      AmErrorAnalyzer(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~AmErrorAnalyzer(); 

      /**
      * Read parameters from archive.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
      
      /**
      * Setup before loop. Opens an output file, if any.
      */
      virtual void setup();

      /**
      * Compute a sampled value and update the accumulator.
      *
      * \param iStep  MD time step index
      */
      virtual void sample(long iStep);

      /**
      * Write final results to file after a simulation.
      */
      virtual void output();
      
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;
      

   protected:
   
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;
      
      /**
      * Write results of statistical analysis to files.
      *
      * \pre hasAccumulator == true
      * \param outputFileName base output file name for analyzer
      */
      void outputAccumulators(std::string outputFileName);      
            
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /**
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();
      
      // Output file stream
      std::ofstream outputFile_;
      
      // Output filename
      std::string filename_;
      
   private:
      
      /// Array of AM prejection step error Average objects for each iteration
      DArray<Average> projectionRatioAccumulators_;
      
      /// Array of AM predicted step error Average objects for each iteration
      DArray<Average> predictRatioAccumulators_;
      
      /// Array of AM mixing step error Average objects for each iteration
      DArray<Average> mixingRatioAccumulators_;
      
      DArray<int> mixingStepCounter_;
      
      /// Average object for itr0 error;
      Average errorItr0Accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      /// Number of maximum iteration.
      int nValue_;
      
      /// Does this processor have accumulators ?
      bool hasAccumulators_;
      
      /// Pointer to parent Simulator.
      Simulator<D>* simulatorPtr_;
      
      /// Pointer to the parent system.
      System<D>* systemPtr_; 

   };
   
   // Inline functions
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& AmErrorAnalyzer<D>::simulator()
   {  return *simulatorPtr_; }
   
   // Get the parent system.
   template <int D>
   inline System<D>& AmErrorAnalyzer<D>::system()
   {  return *systemPtr_; }
   
   #ifndef RPC_AM_ERROR_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class AmErrorAnalyzer<1>;
   extern template class AmErrorAnalyzer<2>;
   extern template class AmErrorAnalyzer<3>;
   #endif
   

}
}
#endif 

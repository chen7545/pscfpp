#ifndef RPC_PERTURBATION_DERIVATIVE_H
#define RPC_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/
#include "Analyzer.h"                             // base class
#include <util/accumulators/Average.h>            // member template

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Compute hamiltonian autocorrelation.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class PerturbationDerivative : public Analyzer<D>
   {

   public:
   
      /**
      * Constructor.
      */
      PerturbationDerivative(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~PerturbationDerivative()
      {} 
      
      /** 
      * Setup before beginning of loop.
      */
      void setup();
   
      /**
      * Sample hamiltonian and add to accumulator.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      void output();
      
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      void readParameters(std::istream& in);
      
      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;     
     
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_; 
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;
      
   protected:      
           
      /**
      * Add current value to accumulator, output if needed.
      *
      * \pre hasAccumulator_== true
      * \param iStep simulation step counter
      */
      void updateAccumulator(long iStep);
      
      
      /**
      * Write results of statistical analysis to files.
      *
      * \pre hasAccumulator_ == true
      * \param outputFileName base output file name for analyzer
      */
      void outputAccumulator(std::string outputFileName);
      
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /** 
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();
 
   private: 
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Output filename
      std::string filename_;
      
      // Output interval 
      long outputInterval_;

      /// Statistical accumulator.
      Average accumulator_;
 
      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      /// Total free energy
      double f_;
      
      /// Has readParam been called?
      bool isInitialized_;
      
      /// Does this processor have accumulators?
      bool hasAccumulator_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& PerturbationDerivative<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& PerturbationDerivative<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_PERTURBATION_DERIVATIVE_TPP
   // Suppress implicit instantiation
   extern template class PerturbationDerivative<1>;
   extern template class PerturbationDerivative<2>;
   extern template class PerturbationDerivative<3>;
   #endif

}
}
#endif 

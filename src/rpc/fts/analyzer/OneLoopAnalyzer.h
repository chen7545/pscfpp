#ifndef RPC_ONE_LOOP_ANALYZER_H
#define RPC_ONE_LOOP_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                             // base class

#include <util/containers/DArray.h>               // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <prdc/cpu/RField.h>                      // member
#include <map>                                    // member

#include <string>
#include <iostream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * OneLoopAnalyzer evaluates gaussian flucutation contribution to 
   * free energy for AB diblock coplymer
   * 
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class OneLoopAnalyzer : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      OneLoopAnalyzer(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      ~OneLoopAnalyzer(){};

      /**
      * Read parameters from file.
      *
      * Input format:
      *   - string            outputFileName  output file base name
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /** 
      * Setup.
      */
      void setup();
      void sample(long iStep){};
   
      /**
      * Output results to predefined output file.
      */
      void output();
      
      /**
      * Compute two point intramolecular function 
      */
      void computeIntra();
      
      /**
      * Compute fluctuation function
      * F(k) = 4S_{++}(k)/ (S_{++}(k)S_{--}(k) - S_{+-}^2(k))
      * S_{++} = S_{AA} + S_{BB} + 2S_{AB}
      * S_{+-} = S_{AA} - S_{BB}
      * S_{--} = S_{AA} + S_{BB} - 2S_{AB}
      */
      void computeFluctuation();
      
   protected:

      /**
      * Output file stream.
      */
      std::ofstream outputFile_;
      
      /**
      * Output filename
      */
      std::string filename_;
      
      /**
      * Intramolecular correlation for S_AA(k),S_BB(k), SAB(k)
      */
      DArray< DArray<double> > intraCorrelations_;
   
      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;     
      
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_; 
      
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /** 
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readOutputFileName;

   private:
      
      /// One loop approximation of fluctuation function
      RField<D> fluctuation_;
      
      /// Size of fourier space
      int kSize_;
      
      /// Total length of diblock copolymer
      double totalN_;
      
      /// Dimensions of fourier space
      IntVec<D> kMeshDimensions_;
      
      /// Has readParam been called?
      bool isInitialized_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& OneLoopAnalyzer<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& OneLoopAnalyzer<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_ONE_LOOP_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class OneLoopAnalyzer<1>;
   extern template class OneLoopAnalyzer<2>;
   extern template class OneLoopAnalyzer<3>;
   #endif

}
}
#endif

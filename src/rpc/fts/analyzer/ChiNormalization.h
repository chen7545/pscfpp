#ifndef RPC_CHI_NORMALIZATION_H
#define RPC_CHI_NORMALIZATION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"

#include <pscf/math/RealVec.h>
#include <prdc/cpu/RFieldDft.h>   
             
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * 
   * The relative number of intermolecular contacts in the limit of 
   * chi_b cpproaches 0 and N approaches to inf.
   * 
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class ChiNormalization : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      ChiNormalization(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~ChiNormalization()
      {}

      /**
      * Read interval and output file name.
      * 
      * Input format:
      *
      *   - string            outputFileName  output file base name
      * 
      * \param in input parameter file
      */
      void readParameters(std::istream& in);
      
      /**
      * Initialization.
      */
      void setup();

      /**
      * Empty function for ChiNormalization analyzer class.
      * 
      * \param iStep step index
      */
      void sample(long iStep){};
      
      /**
      * Output the linear normalization term (zInf).
      */
      void output();
      
      /**
      * Return and compute the probability function.
      * For discrete model, the probability function refers to two 
      * monomoers of a polymer chain that sparated by i bonds occupy s
      * ame grid.
      *
      * \param i two monomer separated by i bonds
      */
      double computeP(int i);
      
      /**
      * Return and compute chi linear normalization term.
      */
      double computeLinearNormalization();
      
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
      * Linear normalization term.
      */
      double zInf_;
      
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
      
      using ParamComposite::readOptional;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readOutputFileName;

   private:

      /// Has readParam been called?
      bool isInitialized_;
      
      /// The cutoff number of bonds 
      int nBond_;
      
      /// Discrete grid spacing along each direction normalized
      /// by statistical length 
      RealVec<D> dx_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& ChiNormalization<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& ChiNormalization<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_CHI_NORMALIZATION_TPP
   // Suppress implicit instantiation
   extern template class ChiNormalization<1>;
   extern template class ChiNormalization<2>;
   extern template class ChiNormalization<3>;
   #endif

}
}
#endif

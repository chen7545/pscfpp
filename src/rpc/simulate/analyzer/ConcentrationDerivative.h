#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/accumulators/Average.h>           // member
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Compute derivative of free energy with respect to chi_bare.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class ConcentrationDerivative : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      ConcentrationDerivative(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~ConcentrationDerivative()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Clear nSample counter.
      */
      virtual void setup();

      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep step index
      */
      virtual void sample(long iStep);

      /**
      * Close trajectory file after run.
      */
      virtual void output();

      using ParamComposite::read;
      using ParamComposite::setClassName;
      using Analyzer<D>::interval;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::isAtInterval;

   protected:

      /**
      * Add current value to accumulator, output if needed.
      *
      * \pre hasAccumulator_== true
      * \param iStep simulation step counter
      */
      void updateAccumulator(long iStep);
      
      /* Write results of statistical analysis to files.
      *
      * \pre hasAccumulator_ == true
      * \param outputFileName base output file name for analyzer
      */
      void outputAccumulator(std::string outputFileName);
      
      /// Output file stream
      std::ofstream outputFile_;

      /// Output filename
      std::string filename_;
      
      /// Statistical accumulator.
      Average accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Has readParam been called?
      long isInitialized_;

      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;

      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;

   protected:

      /**
      * Return reference to parent system.
      */
      System<D>& system();

      /**
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();


   };

   // Get the parent system.
   template <int D>
   inline System<D>& ConcentrationDerivative<D>::system()
   {  return *systemPtr_; }

   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& ConcentrationDerivative<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_CONCENTRATION_DERIVATIVE_TPP
   // Suppress implicit instantiation
   extern template class ConcentrationDerivative<1>;
   extern template class ConcentrationDerivative<2>;
   extern template class ConcentrationDerivative<3>;
   #endif

}
}
#endif

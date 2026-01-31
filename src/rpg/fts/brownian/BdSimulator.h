#ifndef RPG_BD_SIMULATOR_H
#define RPG_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/fts/simulator/Simulator.h>         // base class
#include <rpg/fts/analyzer/AnalyzerManager.h>    // member
//#include <util/param/Factory.h>                  // member template

// Forward declarations
namespace Util {
   template <class T> class Factory;
}
namespace Pscf {
   namespace Rpg {
      template <int D> class System;
      template <int D> class BdStep;
      template <int D> class TrajectoryReader;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Brownian dynamics simulator for PS-FTS.
   *
   * \see \ref rp_BdSimulator_page (manual page)
   *
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class BdSimulator : public Simulator<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      BdSimulator(System<D>& system);

      /**
      * Destructor.
      */
      ~BdSimulator();

      /**
      * Read parameter file block.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Brownian dynamics (BD) simulation.
      *
      * Perform a field theoretic BD simulation using the partial
      * saddle-point approximation.
      *
      * \param nStep  number of BD steps
      */
      void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function creates an instance of the TrajectoryReader subclass
      * specified by the "classname" argument, and uses it to read and
      * analyze a section of a trajectory file, starting at frame number
      * "min" and ending at frame number "max".
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader subclass to use
      * \param filename  name of the trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Does this BdSimulator have a BdStep object?
      */
      bool hasBdStep() const;

      /**
      * Get the BdStep by reference.
      */
      BdStep<D>& stepper();

      /**
      * Get the AnalyzerManger by reference.
      */
      AnalyzerManager<D>& analyzerManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory< TrajectoryReader<D> >& trajectoryReaderFactory();

      ///@}

      // Inherited public functions

      using Simulator<D>::allocate;
      using Simulator<D>::analyzeChi;

      using Simulator<D>::computeWc;
      using Simulator<D>::computeCc;
      using Simulator<D>::computeDc;
      using Simulator<D>::wc;
      using Simulator<D>::cc;
      using Simulator<D>::dc;
      using Simulator<D>::hasWc;
      using Simulator<D>::hasCc;
      using Simulator<D>::hasDc;

      using Simulator<D>::clearData;
      using Simulator<D>::computeHamiltonian;
      using Simulator<D>::hasHamiltonian;
      using Simulator<D>::hamiltonian;
      using Simulator<D>::idealHamiltonian;
      using Simulator<D>::fieldHamiltonian;

      using Simulator<D>::saveState;
      using Simulator<D>::restoreState;
      using Simulator<D>::clearState;
      using Simulator<D>::clearTimers;

      using Simulator<D>::system;
      using Simulator<D>::random;
      using Simulator<D>::vecRandom;
      using Simulator<D>::hasCompressor;
      using Simulator<D>::compressor;
      using Simulator<D>::hasPerturbation;
      using Simulator<D>::perturbation;
      using Simulator<D>::hasRamp;
      using Simulator<D>::ramp;

   protected:

      // Inherited protected functions

      using ParamComposite::readParamComposite;
      using ParamComposite::readOptional;

      using Simulator<D>::readRandomSeed;
      using Simulator<D>::compressorFactory;
      using Simulator<D>::readCompressor;
      using Simulator<D>::perturbationFactory;
      using Simulator<D>::readPerturbation;
      using Simulator<D>::setPerturbation;
      using Simulator<D>::rampFactory;
      using Simulator<D>::readRamp;
      using Simulator<D>::setRamp;

      // Inherited protected data members

      using Simulator<D>::wc_;
      using Simulator<D>::hasWc_;
      using Simulator<D>::hamiltonian_;
      using Simulator<D>::idealHamiltonian_;
      using Simulator<D>::fieldHamiltonian_;
      using Simulator<D>::hasHamiltonian_;
      using Simulator<D>::iStep_;
      using Simulator<D>::iTotalStep_;
      using Simulator<D>::state_;

   private:

      /**
      * Manager for Analyzer.
      */
      AnalyzerManager<D> analyzerManager_;

      /**
      * Pointer to Brownian dynamics step algorithm.
      */
      BdStep<D>* bdStepPtr_;

      /**
      * Pointer to a BdStep factory.
      */
      Factory< BdStep<D> >* bdStepFactoryPtr_;

      /**
      * Pointer to a trajectory reader factory.
      */
      Factory< TrajectoryReader<D> >* trajectoryReaderFactoryPtr_;

      // Private member function

      /**
      * Called at the beginning of the simulation.
      *
      * \param nStep  number of steps planned for the simulation
      */
      void setup(int nStep);

   };

   // Inline member functions

   // Does this BdSimulator have a BdStep?
   template <int D> inline
   bool BdSimulator<D>::hasBdStep() const
   {  return (bool)bdStepPtr_; }

   // Get the BdStep.
   template <int D> inline
   BdStep<D>& BdSimulator<D>::stepper()
   {
      UTIL_CHECK(bdStepPtr_);
      return *bdStepPtr_;
   }

   // Get the analyzer manager.
   template <int D> inline
   AnalyzerManager<D>& BdSimulator<D>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReader factory.
   template <int D> inline
   Factory<TrajectoryReader<D> >& BdSimulator<D>::trajectoryReaderFactory()
   {
      UTIL_CHECK(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

   // Explicit instantiation declarations
   extern template class BdSimulator<1>;
   extern template class BdSimulator<2>;
   extern template class BdSimulator<3>;

}
}
#endif

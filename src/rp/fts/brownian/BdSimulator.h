#ifndef RP_BD_SIMULATOR_H
#define RP_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>          // base class

// Forward declarations
namespace Util {
   template <class T> class Factory;
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Brownian dynamics simulator for PS-FTS.
   *
   * \see \ref rp_BdSimulator_page (manual page)
   *
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class BdSimulator : public Simulator<D, T>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      BdSimulator(typename T::System& system);

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
      typename T::BdStep& bdStep();

      /**
      * Get the AnalyzerManager by reference.
      */
      typename T::AnalyzerManager& analyzerManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory< typename T::TrajectoryReader >& trajectoryReaderFactory();

      ///@}

      using SimulatorT = Simulator<D,T>;

      // Inherited public functions

      using SimulatorT::analyzeChi;
      using SimulatorT::chiEval;
      using SimulatorT::chiEvecs;
      using SimulatorT::computeWc;
      using SimulatorT::computeCc;
      using SimulatorT::computeDc;
      using SimulatorT::wc;
      using SimulatorT::cc;
      using SimulatorT::dc;
      using SimulatorT::hasWc;
      using SimulatorT::hasCc;
      using SimulatorT::hasDc;

      using SimulatorT::clearData;
      using SimulatorT::computeHamiltonian;
      using SimulatorT::hasHamiltonian;
      using SimulatorT::hamiltonian;
      using SimulatorT::idealHamiltonian;
      using SimulatorT::fieldHamiltonian;
      using SimulatorT::perturbationHamiltonian;

      using SimulatorT::saveState;
      using SimulatorT::restoreState;
      using SimulatorT::clearState;
      using SimulatorT::clearTimers;

      using SimulatorT::system;
      using SimulatorT::random;
      using SimulatorT::hasCompressor;
      using SimulatorT::compressor;
      using SimulatorT::hasPerturbation;
      using SimulatorT::perturbation;
      using SimulatorT::hasRamp;
      using SimulatorT::ramp;

   protected:

      // Inherited protected functions

      using ParamComposite::readParamComposite;
      using ParamComposite::readOptional;

      using SimulatorT::readRandomSeed;
      using SimulatorT::compressorFactory;
      using SimulatorT::readCompressor;
      using SimulatorT::perturbationFactory;
      using SimulatorT::readPerturbation;
      using SimulatorT::setPerturbation;
      using SimulatorT::rampFactory;
      using SimulatorT::readRamp;
      using SimulatorT::setRamp;

      // Inherited protected data members

      using SimulatorT::wc_;
      using SimulatorT::hasWc_;
      using SimulatorT::hamiltonian_;
      using SimulatorT::idealHamiltonian_;
      using SimulatorT::fieldHamiltonian_;
      using SimulatorT::hasHamiltonian_;
      using SimulatorT::iStep_;
      using SimulatorT::iTotalStep_;
      using SimulatorT::state_;

   private:

      /**
      * Manager for Analyzer.
      */
      typename T::AnalyzerManager analyzerManager_;

      /**
      * Pointer to Brownian dynamics step algorithm.
      */
      typename T::BdStep* bdStepPtr_;

      /**
      * Pointer to a BdStep factory.
      */
      Factory< typename T::BdStep >* bdStepFactoryPtr_;

      /**
      * Pointer to a trajectory reader factory.
      */
      Factory< typename T::TrajectoryReader >* trajectoryReaderFactoryPtr_;

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
   template <int D, class T> inline
   bool BdSimulator<D,T>::hasBdStep() const
   {  return (bool)bdStepPtr_; }

   // Get the BdStep.
   template <int D, class T> inline
   typename T::BdStep& BdSimulator<D,T>::bdStep()
   {
      UTIL_CHECK(bdStepPtr_);
      return *bdStepPtr_;
   }

   // Get the analyzer manager.
   template <int D, class T> inline
   typename T::AnalyzerManager& BdSimulator<D,T>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReader factory.
   template <int D, class T> inline
   Factory< typename T::TrajectoryReader >& 
   BdSimulator<D,T>::trajectoryReaderFactory()
   {
      UTIL_CHECK(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

}
}
#endif

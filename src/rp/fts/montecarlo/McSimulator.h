#ifndef RP_MC_SIMULATOR_H
#define RP_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>         // base class

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Monte-Carlo simulation coordinator.
   *
   * \see \ref rp_McSimulator_page (Manual Page)
   *
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class McSimulator : public Simulator<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McSimulator(typename T::System& system);

      /**
      * Destructor.
      */
      ~McSimulator();

      /**
      * Read parameter file block.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic Monte-Carlo simulation using the
      * partial saddle-point approximation.
      *
      * \param nStep  number of attempted Monte-Carlo steps
      */
      void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader subclass
      * specified by the "classname" argument to read a trajectory
      * file.
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader subclass to use
      * \param filename  name of the trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      /**
      * Output timing results
      */
      virtual void outputTimers(std::ostream& out) const;

      /**
      * Clear timers
      */
      virtual void clearTimers();

      /**
      * Does the stored state need to include Cc fields?
      */
      bool needsCc();

      /**
      * Does the stored state need to include Dc fields?
      */
      bool needsDc();

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Get the McMoveManager.
      */
      typename T::McMoveManager& mcMoveManager();

      /**
      * Get the AnalyzerManager.
      */
      typename T::AnalyzerManager& analyzerManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory<TrajectoryReader<D>>& trajectoryReaderFactory();

      /**
      * Have any MC moves been defined?
      *
      * Equivalent to a test for mcMoveManager().size() > 0.
      */
      bool hasMcMoves() const;

      ///@}

      using SimulatorT = Simulator<T,D>;

      // Inherited public functions

      using SimulatorT::allocate;
      using SimulatorT::analyzeChi;
      using SimulatorT::chiEval;
      using SimulatorT::chiEvecs;
      using SimulatorT::computeWc;
      using SimulatorT::computeCc;
      using SimulatorT::computeDc;
      using SimulatorT::hasWc;
      using SimulatorT::wc;
      using SimulatorT::hasCc;
      using SimulatorT::cc;
      using SimulatorT::hasDc;
      using SimulatorT::dc;

      using SimulatorT::clearData;
      using SimulatorT::computeHamiltonian;
      using SimulatorT::hasHamiltonian;
      using SimulatorT::hamiltonian;
      using SimulatorT::idealHamiltonian;
      using SimulatorT::fieldHamiltonian;
      using SimulatorT::perturbationHamiltonian;

      using SimulatorT::system;
      using SimulatorT::random;
      using SimulatorT::vecRandom;
      using SimulatorT::hasCompressor;
      using SimulatorT::compressor;
      using SimulatorT::hasPerturbation;
      using SimulatorT::perturbation;
      using SimulatorT::hasRamp;
      using SimulatorT::ramp;

      using SimulatorT::saveState;
      using SimulatorT::restoreState;
      using SimulatorT::clearState;
      using SimulatorT::outputMdeCounter;

   protected:

      // Inherited protected functions

      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::readParamComposite;
      using ParamComposite::readParamCompositeOptional;

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
      * Manager for Monte Carlo moves.
      */
      typename T::McMoveManager mcMoveManager_;

      /**
      * Manager for analyzers.
      */
      typename T::AnalyzerManager analyzerManager_;

      /**
      * Pointer to a trajectory reader factory.
      */
      Factory< TrajectoryReader<D> >* trajectoryReaderFactoryPtr_;

      // Private member function

      /**
      * Setup before the main loop of a simulate or analyze command.
      *
      * \param nStep  number of MC steps to attempt
      */
      void setup(int nStep);

   };

   // Get the Monte-Carlo move manager.
   template <int D, class T> inline
   typename T::McMoveManager& McSimulator<D,T>::mcMoveManager()
   {  return mcMoveManager_; }

   // Get the analyzer manager.
   template <int D, class T> inline
   typename T::AnalyzerManager& McSimulator<D,T>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReader factory.
   template <int D, class T> inline
   Factory<TrajectoryReader<D> >& McSimulator<D,T>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

   // Have any MC moves been defined?
   template <int D, class T> inline
   bool McSimulator<D,T>::hasMcMoves() const
   {  return (bool)(mcMoveManager_.size() > 0); }

   // Does the stored state need to include Cc fields?
   template <int D, class T> inline
   bool McSimulator<D,T>::needsCc()
   {  return state_.needsCc; }

   // Does the stored state need to include Dc fields?
   template <int D, class T> inline
   bool McSimulator<D,T>::needsDc()
   {  return state_.needsDc; }

}
}
#endif

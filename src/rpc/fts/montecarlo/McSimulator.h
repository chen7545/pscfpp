#ifndef RPC_MC_SIMULATOR_H
#define RPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/simulator/Simulator.h>         // base class
#include <rpc/fts/montecarlo/McMoveManager.h>    // member
#include <rpc/fts/analyzer/AnalyzerManager.h>    // member

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class McMove;
   template <int D> class TrajectoryReader;

   /**
   * Monte-Carlo simulation coordinator.
   *
   * \see \ref rp_McSimulator_page (Manual Page)
   *
   * \ingroup Rpc_Fts_MonteCarlo_Module
   */
   template <int D>
   class McSimulator : public Simulator<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McSimulator(System<D>& system);

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
      McMoveManager<D>& mcMoveManager();

      /**
      * Get the AnalyzerManager.
      */
      AnalyzerManager<D>& analyzerManager();

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

      // Inherited public functions

      using Simulator<D>::allocate;
      using Simulator<D>::analyzeChi;
      using Simulator<D>::chiEval;
      using Simulator<D>::chiEvecs;
      using Simulator<D>::computeWc;
      using Simulator<D>::computeCc;
      using Simulator<D>::computeDc;
      using Simulator<D>::hasWc;
      using Simulator<D>::wc;
      using Simulator<D>::hasCc;
      using Simulator<D>::cc;
      using Simulator<D>::hasDc;
      using Simulator<D>::dc;

      using Simulator<D>::clearData;
      using Simulator<D>::computeHamiltonian;
      using Simulator<D>::hasHamiltonian;
      using Simulator<D>::hamiltonian;
      using Simulator<D>::idealHamiltonian;
      using Simulator<D>::fieldHamiltonian;
      using Simulator<D>::perturbationHamiltonian;

      using Simulator<D>::system;
      using Simulator<D>::random;
      using Simulator<D>::vecRandom;
      using Simulator<D>::hasCompressor;
      using Simulator<D>::compressor;
      using Simulator<D>::hasPerturbation;
      using Simulator<D>::perturbation;
      using Simulator<D>::hasRamp;
      using Simulator<D>::ramp;

      using Simulator<D>::saveState;
      using Simulator<D>::restoreState;
      using Simulator<D>::clearState;
      using Simulator<D>::outputMdeCounter;

   protected:

      // Inherited protected functions

      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::readParamComposite;
      using ParamComposite::readParamCompositeOptional;

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
      * Manager for Monte Carlo moves.
      */
      McMoveManager<D> mcMoveManager_;

      /**
      * Manager for analyzers.
      */
      AnalyzerManager<D> analyzerManager_;

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
   template <int D> inline
   McMoveManager<D>& McSimulator<D>::mcMoveManager()
   {  return mcMoveManager_; }

   // Get the analyzer manager.
   template <int D> inline
   AnalyzerManager<D>& McSimulator<D>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReader factory.
   template <int D> inline
   Factory<TrajectoryReader<D> >& McSimulator<D>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

   // Have any MC moves been defined?
   template <int D> inline
   bool McSimulator<D>::hasMcMoves() const
   {  return (bool)(mcMoveManager_.size() > 0); }

   // Does the stored state need to include Cc fields?
   template <int D> inline
   bool McSimulator<D>::needsCc()
   {  return state_.needsCc; }

   // Does the stored state need to include Dc fields?
   template <int D> inline
   bool McSimulator<D>::needsDc()
   {  return state_.needsDc; }

   // Explicit instantiation declarations
   extern template class McSimulator<1>;
   extern template class McSimulator<2>;
   extern template class McSimulator<3>;

}
}
#endif

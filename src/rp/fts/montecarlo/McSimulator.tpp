#ifndef RP_MC_SIMULATOR_TPP
#define RP_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <util/param/Factory.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/misc/Timer.h>

#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   McSimulator<D,T>::McSimulator(SystemT& system, McSimulatorT& mcSimulator)
    : SimulatorT(system),
      mcMoveManager_(mcSimulator, system),
      analyzerManager_(mcSimulator, system),
      trajectoryReaderFactoryPtr_(nullptr)
   {
      ParamComposite::setClassName("McSimulator");
      trajectoryReaderFactoryPtr_
             = new typename T::TrajectoryReaderFactory(system);
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   McSimulator<D,T>::~McSimulator()
   {
      if (trajectoryReaderFactoryPtr_) {
         delete trajectoryReaderFactoryPtr_;
      }
   }

   /*
   * Read parameter file block.
   */
   template <int D, class T>
   void McSimulator<D,T>::readParameters(std::istream &in)
   {
      // Optionally read random seed. Initialize random number generators
      SimulatorT::readRandomSeed(in);

      // Optionally read McMoveManager block
      ParamComposite::readParamCompositeOptional(in, mcMoveManager_);

      // Optionally read Compressor block
      bool isEnd = false;
      SimulatorT::readCompressor(in, isEnd);
      if (hasMcMoves()) {
         UTIL_CHECK(hasCompressor());
      }

      // A Compressor is required if MC moves are declared.
      // A Ramp is allowed only if MC moves are declared.

      // Optionally read Perturbation and/or Ramp blocks
      SimulatorT::readPerturbation(in, isEnd);
      if (hasMcMoves()) {
         SimulatorT::readRamp(in, isEnd);
      }

      // Optionally read AnalyzerManager block
      AnalyzerT::baseInterval = 0; // default value
      ParamComposite::readParamCompositeOptional(in, analyzerManager_);

      // Figure out what needs to be saved in stored state
      state_.needsCc = false;
      state_.needsDc = false;
      state_.needsHamiltonian = true;
      if (mcMoveManager_.needsCc()){
         state_.needsCc = true;
      }
      if (mcMoveManager_.needsDc()){
         state_.needsDc = true;
      }

      // Allocate memory for SimulatorT base class
      SimulatorT::allocate();
   }

   /*
   * Initialize just prior to a run.
   */
   template <int D, class T>
   void McSimulator<D,T>::setup(int nStep)
   {
      UTIL_CHECK(SimulatorT::system().w().hasData());

      // Eigenanalysis of the projected chi matrix.
      SimulatorT::analyzeChi();

      if (hasPerturbation()) {
         SimulatorT::perturbation().setup();
      }

      if (hasRamp()) {
         SimulatorT::ramp().setup(nStep);
      }

      // Solve the MDE and compute c-fields for the initial state
      SimulatorT::system().compute();

      // Compress the initial state (adjust pressure-like field)
      if (hasCompressor()) {
         SimulatorT::compressor().compress();
         SimulatorT::compressor().clearTimers();
      }

      // Compute field components and Hamiltonian
      SimulatorT::computeWc();
      if (state_.needsCc || state_.needsDc) {
         SimulatorT::computeCc();
      }
      if (state_.needsDc) {
         SimulatorT::computeDc();
      }
      SimulatorT::computeHamiltonian();

      // Setup MC moves (if any)
      if (hasMcMoves()) {
         mcMoveManager_.setup();
      }

      // Setup analyzers (if any)
      if (analyzerManager_.size() > 0){
         analyzerManager_.setup();
      }

   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D, class T>
   void McSimulator<D,T>::simulate(int nStep)
   {
      UTIL_CHECK(hasMcMoves());
      UTIL_CHECK(hasCompressor());

      // Initial setup
      setup(nStep);
      iStep_ = 0;
      if (hasRamp()) {
         SimulatorT::ramp().setParameters(iStep_);
      }
      int analyzerBaseInterval = AnalyzerT::baseInterval;
      Log::file() << std::endl;

      // Start timers
      Timer timer;
      Timer analyzerTimer;
      timer.start();

      // Analyze initial step
      analyzerTimer.start();
      analyzerManager_.sample(iStep_);
      analyzerTimer.stop();

      // Main Monte Carlo loop
      for (iTotalStep_ = 0; iTotalStep_ < nStep; ++iTotalStep_) {

         // Choose and attempt an McMove
         bool converged;
         converged = mcMoveManager_.chooseMove().move();

         if (converged){
            iStep_++;

            if (hasRamp()) {
               SimulatorT::ramp().setParameters(iStep_);
            }

            // Analysis (if any)
            analyzerTimer.start();
            if (analyzerBaseInterval != 0) {
               if (iStep_ % analyzerBaseInterval == 0) {
                  if (analyzerManager_.size() > 0) {
                     analyzerManager_.sample(iStep_);
                  }
               }
            }
            analyzerTimer.stop();

         } else{
            Log::file() << "Step: "<< iTotalStep_
                        << " failed to converge" << "\n";
         }
      }
      timer.stop();
      double time = timer.time();
      double analyzerTime = analyzerTimer.time();

      // Output move statistics and final analysis results
      mcMoveManager_.output();
      if (analyzerBaseInterval != 0){
         analyzerManager_.output();
      }

      // Final output from ramp (if any)
      if (SimulatorT::hasRamp()){
         SimulatorT::ramp().output();
      }

      // Output times for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep               " << nStep << std::endl;
      if (iStep_ != nStep){
         Log::file() << "nFail Step          " 
                     << (nStep - iStep_) << std::endl;
      }
      Log::file() << "Total run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep        " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << "Analyzer run time   " << analyzerTime
                  << " sec" << std::endl;
      Log::file() << std::endl;

      // Output number of times MDE has been solved for the simulation run
      outputMdeCounter(Log::file());

      // Print McMove acceptance statistics
      long attempt;
      long accept;
      long fail;
      using namespace std;
      Log::file() << "Move Statistics:" << endl << endl;
      Log::file() << setw(20) << left <<  "Move Name"
           << setw(10) << right << "Attempted"
           << setw(10) << right << "Accepted"
           << setw(13) << right << "AcceptRate"
           << setw(10) << right << "Failed"
           << setw(13) << right << "FailRate"
           << endl;
      int nMove = mcMoveManager_.size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = mcMoveManager_[iMove].nAttempt();
         accept  = mcMoveManager_[iMove].nAccept();
         fail  = mcMoveManager_[iMove].nFail();
         Log::file() << setw(20) << left
              << mcMoveManager_[iMove].className()
              << setw(10) << right << attempt
              << setw(10) << accept
              << setw(13) << fixed << setprecision(5)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << setw(10) << fail
              << setw(13) << fixed << setprecision(5)
              << ( attempt == 0 ? 0.0 : double(fail)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

   }

   /*
   * Open, read and analyze a trajectory file
   */
   template <int D, class T>
   void McSimulator<D,T>::analyze(int min, int max,
                                  std::string classname,
                                  std::string filename)
   {
      // Preconditions
      UTIL_CHECK(min >= 0);
      UTIL_CHECK(max >= min);
      UTIL_CHECK(AnalyzerT::baseInterval > 0);
      UTIL_CHECK(analyzerManager_.size() > 0);

      // Construct TrajectoryReader
      typename T::TrajectoryReader* trajectoryReaderPtr;
      trajectoryReaderPtr = trajectoryReaderFactory().factory(classname);
      if (!trajectoryReaderPtr) {
         std::string message;
         message = "Invalid TrajectoryReader class name " + classname;
         UTIL_THROW(message.c_str());
      }

      // Open trajectory file
      trajectoryReaderPtr->open(filename);
      trajectoryReaderPtr->readHeader();

      // Main loop over trajectory frames
      Timer timer;
      timer.start();
      bool hasFrame = trajectoryReaderPtr->readFrame();

      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         if (hasFrame) {
            SimulatorT::clearData();

            // Initialize analyzers
            if (iStep_ == min) {
               setup(iStep_);
            }

            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager_.sample(iStep_);
            }
         }

         hasFrame = trajectoryReaderPtr->readFrame();
      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;
      int nFrames = iStep_ - min;
      trajectoryReaderPtr->close();
      delete trajectoryReaderPtr;

      // Output results of all analyzers to output files
      analyzerManager_.output();

      // Output number of frames and times
      Log::file() << std::endl;
      Log::file() << "# of frames   " << nFrames << std::endl;
      Log::file() << "run time      " << timer.time()
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames)
                  << "  sec" << std::endl;
      Log::file() << std::endl;

   }

   /*
   * Output McMoveManager timer results.
   */
   template <int D, class T>
   void McSimulator<D,T>::outputTimers(std::ostream& out) const
   {
      SimulatorT::compressor().outputTimers(out);
      out << "\n";
      out << "MC move time contributions:\n";
      mcMoveManager_.outputTimers(out);
   }

   /*
   * Clear all McMoveManager timers.
   */
   template <int D, class T>
   void McSimulator<D,T>::clearTimers()
   {  mcMoveManager_.clearTimers(); }

}
}
#endif

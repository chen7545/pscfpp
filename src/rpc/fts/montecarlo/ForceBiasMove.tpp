#ifndef RPC_FORCE_BIAS_MOVE_TPP
#define RPC_FORCE_BIAS_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMove.h"
#include "McMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/system/System.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>



namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D>::ForceBiasMove(McSimulator<D>& simulator)
    : McMove<D>(simulator),
      w_(),
      dwc_(),
      mobility_(0.0)
   {  ParamComposite::setClassName("ForceBiasMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   ForceBiasMove<D>::~ForceBiasMove()
   {}

   /*
   * Read body of parameter file block and allocate memory.
   */
   template <int D>
   void ForceBiasMove<D>::readParameters(std::istream &in)
   {
      McMoveT::readProbability(in);
      ParamComposite::read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dc_.allocate(nMonomer-1);
      dwc_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         dc_[i].allocate(meshDimensions);
         dwc_[i].allocate(meshDimensions);
      }
      biasField_.allocate(meshDimensions);
      eta_.allocate(meshDimensions);
   }

   /*
   * Setup before entering main simulation loop.
   */
   template <int D>
   void ForceBiasMove<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dc_.capacity() == nMonomer-1);
      UTIL_CHECK(dwc_.capacity() == nMonomer-1);
      for (int i = 0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dc_[i].capacity() == meshSize);
         UTIL_CHECK(dwc_[i].capacity() == meshSize);
      }

      McMove<D>::setup();
   }

   /*
   * Attempt and accept or reject MC move
   */
   template <int D>
   bool ForceBiasMove<D>::move()
   {
      McMoveT::totalTimer_.start();
      McMoveT::incrementNAttempt();

      // Preconditions
      UTIL_CHECK(simulator().hasWc());
      UTIL_CHECK(simulator().hasDc());
      UTIL_CHECK(simulator().hasHamiltonian());

      // Array sizes and index variables
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;

      // Get current Hamiltonian
      double oldHamiltonian = simulator().hamiltonian();

      // Save current state
      simulator().saveState();

      // Clear eigen-components of the fields and Hamiltonian
      simulator().clearData();

      McMoveT::attemptMoveTimer_.start();

      // Copy current W fields from parent system into wc_
      for (i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Copy derivative fields into dc_
      for (i = 0; i < nMonomer - 1; ++i) {
         VecOp::eqV(dc_[i], simulator().dc(i));
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double vNode = vSystem/double(meshSize);
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_/vNode);
      const double stddev = 1.0;
      const double mean = 0.0;

      // Modify local variables dwc_ and w_
      // Loop over eigenvectors of projected chi matrix
      for (j = 0; j < nMonomer - 1; ++j) {

         // Generate a vector of normal distributed random numbers
         McMoveT::vecRandom().normal(eta_, stddev, mean);

         // Compute vector dwc_[j] of field component changes
         VecOp::addVcVc(dwc_[j], dc_[j], a, eta_, b);

         // Loop over monomer types to add to w_
         for (i = 0; i < nMonomer; ++i) {
            double evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(w_[i], dwc_[j], evec);
         }

      }

      // Set modified fields in parent system
      system().w().setRGrid(w_);
      simulator().clearData();

      McMoveT::attemptMoveTimer_.stop();

      // Call compressor
      McMoveT::compressorTimer_.start();
      int compress = simulator().compressor().compress();
      McMoveT::compressorTimer_.stop();

      bool isConverged = false;
      if (compress != 0){
         McMoveT::incrementNFail();
         simulator().restoreState();
      } else {
         isConverged = true;

         // Compute eigenvector components of current fields
         McMoveT::componentTimer_.start();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();
         McMoveT::componentTimer_.stop();

         // Evaluate new Hamiltonian
         McMoveT::hamiltonianTimer_.start();
         simulator().computeHamiltonian();
         double newHamiltonian = simulator().hamiltonian();
         double dH = newHamiltonian - oldHamiltonian;

         // Compute force bias
         double bias = 0.0;
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D> const & di = dc_[j];
            RField<D> const & df = simulator().dc(j);
            computeForceBias(biasField_, di, df, dwc_[j], mobility_);
            bias += Reduce::sum(biasField_);
         }
         bias *= vNode;
         McMoveT::hamiltonianTimer_.stop();

         // Accept or reject move
         bool accept = false;
         McMoveT::decisionTimer_.start();
         double weight = exp(bias - dH);
         accept = McMoveT::random().metropolis(weight);
         if (accept) {
            McMoveT::incrementNAccept();
            simulator().clearState();
         } else {
            simulator().restoreState();
         }
         McMoveT::decisionTimer_.stop();

      }

      McMoveT::totalTimer_.stop();
      return isConverged;
   }

   /*
   * Output data to log file (do-nothing implementation).
   */
   template <int D>
   void ForceBiasMove<D>::output()
   {}

   template<int D>
   void ForceBiasMove<D>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << "ForceBiasMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

   // Private member function

   /*
   * Compute force bias field for use in Metropolis acceptance test.
   */
   template<int D>
   void ForceBiasMove<D>::computeForceBias(
                               RField<D>& result, 
                               RField<D> const & di, 
                               RField<D> const & df, 
                               RField<D> const & dwc, 
                               double mobility)
   {
      const int n = system().domain().mesh().size();
      UTIL_CHECK(result.capacity() == n);
      UTIL_CHECK(di.capacity() == n);
      UTIL_CHECK(df.capacity() == n);
      UTIL_CHECK(dwc.capacity() == n);

      double dp, dm;
      for (int k = 0; k < n; ++k) {
         dp = 0.5*(di[k] + df[k]);
         dm = 0.5*(di[k] - df[k]);
         result[k] = dp*( dwc[k] + mobility*dm );
      }
   }

}
}
#endif

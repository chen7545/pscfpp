#ifndef RPC_REAL_MOVE_TPP
#define RPC_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealMove<D>::RealMove(McSimulator<D>& simulator)
    : McMove<D>(simulator),
      dwc_(),
      sigma_(0.0),
      isAllocated_(false)
   {  ParamComposite::setClassName("RealMove"); }

   /*
   * Destructor.
   */
   template <int D>
   RealMove<D>::~RealMove()
   {}

   /*
   * Read body of parameter file block.
   */
   template <int D>
   void RealMove<D>::readParameters(std::istream &in)
   {
      McMove<D>::readProbability(in);

      // Standard deviation of field change
      ParamComposite::read(in, "sigma", sigma_);
   }

   /*
   * Setup before simulation loop.
   */
   template <int D>
   void RealMove<D>::setup()
   {
      McMove<D>::setup();
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & meshDimensions = system().domain().mesh().dimensions();
      if (!isAllocated_){
         dwc_.allocate(meshDimensions);
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(meshDimensions);
         }
         isAllocated_ = true;
      }
   }

   /*
   * Attempt unconstrained move
   */
   template <int D>
   void RealMove<D>::attemptMove()
   {
      // Copy current fields to w_
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         w_[i] = system().w().rgrid(i);
      }

      // Loop over composition eigenvectors of projected chi matrix
      double evec, mean;
      mean = 0.0;
      for (int j = 0; j < nMonomer - 1; j++){

         // Generate random field changes
         vecRandom().normal(dwc_, sigma_, mean);

         // Add changes to w_ field components
         for (int i = 0; i < nMonomer; ++i) {
            evec = simulator().chiEvecs(j, i);
            VecOp::addEqVc(w_[i], dwc_, evec);
         }
      }

      // Update w-fields in parent system
      system().w().setRGrid(w_);
   }

   template<int D>
   void RealMove<D>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << "RealMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif

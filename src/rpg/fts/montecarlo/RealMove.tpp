#ifndef RPG_REAL_MOVE_TPP
#define RPG_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h"
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealMove<D>::RealMove(McSimulator<D>& simulator)
    : McMove<D>(simulator),
      w_(),
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

      // Standard deviation of field changes
      ParamComposite::read(in, "sigma", sigma_);
   }

   /*
   * Setup before simulation loop.
   */
   template <int D>
   void RealMove<D>::setup()
   {
      McMove<D>::setup();

      if (!isAllocated_){
         const int nMonomer = system().mixture().nMonomer();
         IntVec<D> meshDimensions = system().domain().mesh().dimensions();
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(meshDimensions);
         }
         dwc_.allocate(meshDimensions);
         isAllocated_ = true;
      }
   }

   /*
   * Attempt unconstrained move
   */
   template <int D>
   void RealMove<D>::attemptMove()
   {
      // Copy current W fields to w_
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Loop over composition eigenvectors of projected chi matrix
      double evec, mean;
      mean = 0.0;
      for (int j = 0; j < nMonomer - 1; ++j) {

         // Generate random field changes
         vecRandom().normal(dwc_, sigma_, mean);

         // Add changes to w_ field components
         for (int i = 0; i < nMonomer; ++i) {
            evec = simulator().chiEvecs(j, i);
            VecOp::addEqVc(w_[i], dwc_, evec);
         }
      }

      // Update w-fields of parent system
      system().w().setRGrid(w_);
   }

   /*
   * Output time contributions.
   */
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

#ifndef RPC_LR_AM_COMPRESSOR_TPP
#define RPC_LR_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <prdc/cpu/FFT.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>


namespace Pscf {
namespace Rpc{

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
    : Compressor<D>(system),
      intra_(system),
      isIntraCalculated_(false),
      isAllocated_(false)
   {  ParamComposite::setClassName("LrAmCompressor"); }

   /*
   * Destructor.
   */
   template <int D>
   LrAmCompressor<D>::~LrAmCompressor()
   {}

   /*
   * Read parameters from file.
   */
   template <int D>
   void LrAmCompressor<D>::readParameters(std::istream& in)
   {
      // Default values
      AmTmpl::maxItr_ = 100;
      AmTmpl::verbose_ = 0;
      AmTmpl::errorType_ = "rms";

      // Call base class read methods
      AmTmpl::readParameters(in);
      AmTmpl::readErrorType(in);
   }

   /*
   * Initialize just before entry to iterative loop.
   */
   template <int D>
   void LrAmCompressor<D>::setup(bool isContinuation)
   {
      
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DArray<double> >::setup(isContinuation);
      
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
     
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate memory required by compressor, if not done previously
      if (!isAllocated_){
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         intraCorrelationK_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(dimensions);
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }
      
      // Compute intraCorrelation
      if (!isIntraCalculated_){
         intra_.computeIntraCorrelations(intraCorrelationK_);
         isIntraCalculated_ = true;
      }
      
      // Store initial values of monomer chemical potential fields
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0_[i], system().w().rgrid(i));
      }
      
   }

   /*
   * Main function - iteratively adjust the pressure field.
   */
   template <int D>
   int LrAmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DArray<double> >::solve();
      return solve;
   }

   /*
   * Output timer results.
   */
   template<int D>
   void LrAmCompressor<D>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "LrAmCompressor time contributions:\n";
      AmTmpl::outputTimers(out);
   }

   /*
   * Clear timers and the MDE counter
   */
   template<int D>
   void LrAmCompressor<D>::clearTimers()
   {
      AmTmpl::clearTimers();
      Compressor<D>::mdeCounter_ = 0;
   }

   // Private virtual AM algorithm operations

   /*
   * Correction step (second step of Anderson mixing)
   *
   * This LrAM algorithm uses a quasi-Newton correction step with an
   * approximate Jacobian given by the Jacobian in a homogeneous state.
   */
   template <int D>
   void
   LrAmCompressor<D>::addCorrection(DArray<double>& fieldTrial,
                                    DArray<double> const & resTrial)
   {
      // Copy resTrial to RField<D> resid_ 
      // Allows use of FFT functions that take an RField container
      VecOp::eqV(resid_, resTrial);

      // Convert residual to Fourier space
      system().domain().fft().forwardTransform(resid_, residK_);

      // Combine with linear response factor to obtain update step
      const double vMonomer = system().mixture().vMonomer();
      VecOp::divEqVc(residK_, intraCorrelationK_, vMonomer);

      // Convert update back to real space (destroys residK_)
      system().domain().fft().inverseTransformUnsafe(residK_, resid_);

      // Add update to obtain new fieldTrial
      VecOp::addEqV(fieldTrial, resid_);
   }

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D>
   int LrAmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Does the associated system have initialized w fields?
   */
   template <int D>
   bool LrAmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field variable from the system.
   *
   * The field variable is the change in the Lagrange multiplier field.
   * relative to that used in initial array of monomer fields, w0_.
   */
   template <int D>
   void LrAmCompressor<D>::getCurrent(DArray<double>& curr)
   {
      /*
      * The field that we are adjusting is the Langrange multiplier 
      * field.  The current value is the difference between w and w0_ 
      * for the first monomer type, but any monomer type would give
      * the same answer.
      */
      VecOp::subVV(curr, system().w().rgrid(0), w0_[0]);
   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void LrAmCompressor<D>::evaluate()
   {
      system().compute();
      ++(Compressor<D>::mdeCounter_);
   }

   /*
   * Compute the residual vector for the current system state
   */
   template <int D>
   void LrAmCompressor<D>::getResidual(DArray<double>& resid)
   {

      // Initialize residual to -1.0
      VecOp::eqS(resid, -1.0);

      // Add c fields to get SCF residual vector
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::addEqV(resid, system().c().rgrid(i));
      }
   }

   /*
   * Update the w field values stored in the system 
   */
   template <int D>
   void LrAmCompressor<D>::update(DArray<double>& newGuess)
   {
      // New field is w0 + newGuess for the pressure field
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(wFieldTmp_[i], w0_[i], newGuess);
      }

      // Update system r-grid fields
      system().w().setRGrid(wFieldTmp_);
   }

   /*
   * Output results to log file (do-nothing implementation).
   */
   template<int D>
   void LrAmCompressor<D>::outputToLog()
   {}

}
}
#endif

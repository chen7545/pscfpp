#ifndef RPG_LR_AM_COMPRESSOR_TPP
#define RPG_LR_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/cuda/FFT.h>
#include <prdc/cuda/resources.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
    : intra_(system),
      isIntraCalculated_(false),
      isAllocated_(false)
   {
      CompressorT::setSystem(system);  
      ParamComposite::setClassName("LrAmCompressor"); 
   }

   /*
   * Destructor.
   */
   template <int D>
   LrAmCompressor<D>::~LrAmCompressor()
   {}

   /*
   * Read body of parameter file block and initialize.
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

      // Allocate memory required by AM algorithm, if not done earlier
      AmTmpl::setup(isContinuation);

      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate memory required by compressor, if not done earlier
      if (!isAllocated_) {
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

      // Compute intraCorrelationK_
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
      int solve = AmTmpl::solve();
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
   * Clear timers and the MDE counter.
   */
   template<int D>
   void LrAmCompressor<D>::clearTimers()
   {
      AmTmpl::clearTimers();
      CompressorT::mdeCounter_ = 0;
   }

   // Private virtual AM algorithm operation functions

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D>
   int LrAmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Does the system have an initial field guess?
   */
   template <int D>
   bool LrAmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field variable from the system.
   *
   * The field variable is the change in the Lagrange multiplier field
   * relative to that used in the initial array of monomer fields, w0_.
   */
   template <int D>
   void LrAmCompressor<D>::getCurrent(VectorT& curr)
   {
      /*
      * The field that we are adjusting is the Langrange multiplier
      * field. The current value is the difference between w and w0_
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
      ++(CompressorT::mdeCounter_);
   }

   /*
   * Compute the residual vector for the current system state.
   */
   template <int D>
   void LrAmCompressor<D>::getResidual(VectorT& resid)
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
   * Correction step (second step of Anderson mixing)
   *
   * This LrAm algorithm uses a quasi-Newton correction step with an
   * approximate Jacobian given by the Jacobian in a homogeneous state.
   */
   template <int D>
   void
   LrAmCompressor<D>::addCorrection(VectorT& fieldTrial,
                                    VectorT const & resTrial)
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

      // Add update resid_ to obtain corrected fieldTrial_
      VecOp::addEqV(fieldTrial, resid_);
   }

   /*
   * Update the w field values stored in the system
   */
   template <int D>
   void LrAmCompressor<D>::update(VectorT& newGuess)
   {
      // New field is w0_ + newGuess for the pressure field
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
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

   // Private virtual functions for vector math

   /*
   * Vector assignment, a = b.
   */
   template <int D>
   void LrAmCompressor<D>::setEqual(VectorT& a,
                                    VectorT const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b);
   }

   /*
   * Compute and return inner product of two vectors.
   */
   template <int D>
   double LrAmCompressor<D>::dotProduct(VectorT const & a,
                                        VectorT const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   /*
   * Compute and return maximum element of a vector.
   */
   template <int D>
   double LrAmCompressor<D>::maxAbs(VectorT const & a)
   {  return Reduce::maxAbs(a); }

   /*
   * Compute the vector difference a = b - c
   */
   template <int D>
   void LrAmCompressor<D>::subVV(VectorT& a,
                                 VectorT const & b,
                                 VectorT const & c)
   {  VecOp::subVV(a, b, c); }

   /*
   * Compute a += b*c for vectors a and b, scalar c
   */
   template <int D>
   void LrAmCompressor<D>::addEqVc(VectorT& a,
                                   VectorT const & b,
                                   double c)
   {  VecOp::addEqVc(a, b, c); }

}
}
#endif

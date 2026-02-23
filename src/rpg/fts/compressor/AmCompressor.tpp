#ifndef RPG_AM_COMPRESSOR_TPP
#define RPG_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/cuda/RField.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
    : Compressor<D>(system),
      isAllocated_(false)
   {  ParamComposite::setClassName("AmCompressor"); }

   /*
   * Destructor.
   */
   template <int D>
   AmCompressor<D>::~AmCompressor()
   {}

   /*
   * Read parameters from file.
   */
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      // Default values
      AmTmpl::maxItr_ = 100;
      AmTmpl::verbose_ = 0;
      AmTmpl::errorType_ = "rms";
      bool useLambdaRamp = false;

      AmTmpl::readParameters(in);
      AmTmpl::readErrorType(in);
      AmTmpl::readMixingParameters(in, useLambdaRamp);
   }

   /*
   * Initialize just before entry to iterative loop.
   */
   template <int D>
   void AmCompressor<D>::setup(bool isContinuation)
   {
      // Allocate memory required by AM algorithm if not done earlier.
      AmTmpl::setup(isContinuation);

      const int nMonomer = system().mixture().nMonomer();
      const IntVec<D> dimensions = system().domain().mesh().dimensions();

      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_) {
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(dimensions);
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }

      // Store value of initial guess chemical potential fields
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0_[i], system().w().rgrid(i));
      }
   }

   /*
   * Main function - identify partial saddle-point state.
   */
   template <int D>
   int AmCompressor<D>::compress()
   {
      int solve = AmTmpl::solve();
      return solve;
   }

   /*
   * Output timer information, if requested.
   */
   template<int D>
   void AmCompressor<D>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "Compressor time contributions:\n";
      AmTmpl::outputTimers(out);
   }

   /*
   * Clear timers and MDE counter.
   */
   template<int D>
   void AmCompressor<D>::clearTimers()
   {
      AmTmpl::clearTimers();
      Compressor<D>::mdeCounter_ = 0;
   }

   // Private virtual functions that interact with the parent System

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D>
   int AmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Does the system have an initial field guess?
   */
   template <int D>
   bool AmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field from the system.
   */
   template <int D>
   void AmCompressor<D>::getCurrent(VectorT& curr)
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
   void AmCompressor<D>::evaluate()
   {
      system().compute();
      ++(Compressor<D>::mdeCounter_);
   }

   /*
   * Compute the residual vector for the current system state.
   */
   template <int D>
   void AmCompressor<D>::getResidual(VectorT& resid)
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
   * Update the current system field coordinates.
   */
   template <int D>
   void AmCompressor<D>::update(VectorT& newGuess)
   {
      // New field is w0_ + newGuess for the pressure field
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::addVV(wFieldTmp_[i], w0_[i], newGuess);
      }

      // Set system r-grid fields
      system().w().setRGrid(wFieldTmp_);
   }

   /*
   * Do-nothing output function.
   */
   template<int D>
   void AmCompressor<D>::outputToLog()
   {}

   // Private virtual functions for vector math

   /*
   * Assign one array to another.
   */
   template <int D>
   void AmCompressor<D>::setEqual(VectorT& a,
                                  VectorT const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b);
   }

   /*
   * Compute and return inner product of two vectors.
   */
   template <int D>
   double AmCompressor<D>::dotProduct(VectorT const & a,
                                      VectorT const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   /*
   * Compute and return maximum absolute value element of a vector.
   */
   template <int D>
   double AmCompressor<D>::maxAbs(VectorT const & a)
   {  return Reduce::maxAbs(a); }

   /*
   * Compute the vector difference a = b - c
   */
   template <int D>
   void AmCompressor<D>::subVV(VectorT& a,
                               VectorT const & b,
                               VectorT const & c)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      VecOp::subVV(a, b, c);
   }

   /*
   * Composite a += b*c for vectors a and b, scalar c
   */
   template <int D>
   void AmCompressor<D>::addEqVc(VectorT& a,
                                 VectorT const & b,
                                 double c)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      VecOp::addEqVc(a, b, c);
   }

}
}
#endif

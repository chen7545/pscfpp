#ifndef RP_AM_COMPRESSOR_TPP
#define RP_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpg/field/Domain.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D, class T, class V>
   AmCompressor<D,T,V>::AmCompressor(typename T::System& system)
    : CompressorT(system),
      isAllocated_(false)
   {  ParamComposite::setClassName("AmCompressor"); }

   /*
   * Destructor.
   */
   template <int D, class T, class V>
   AmCompressor<D,T,V>::~AmCompressor()
   {}

   /*
   * Read parameters from file.
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::readParameters(std::istream& in)
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
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::setup(bool isContinuation)
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
   * Main compressor function.
   */
   template <int D, class T, class V>
   int AmCompressor<D,T,V>::compress()
   {
      int solve = AmTmpl::solve();
      return solve;
   }

   /*
   * Output timer information, if requested.
   */
   template<int D>
   void AmCompressor<D,T,V>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "Compressor times contributions:\n";
      AmTmpl::outputTimers(out);
   }

   /*
   * Clear timers and MDE counter.
   */
   template<int D>
   void AmCompressor<D,T,V>::clearTimers()
   {
      AmTmpl::clearTimers();
      CompressorT::mdeCounter_ = 0;
   }

   // Private virtual functions that interact with the parent System

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D, class T, class V>
   int AmCompressor<D,T,V>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Does the system have an initial field guess?
   */
   template <int D, class T, class V>
   bool AmCompressor<D,T,V>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field from the system.
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::getCurrent(V& curr)
   {
      /*
      * The field that we are adjusting is the Langrange multiplier
      * field.  The current value is the difference between w and w0_
      * for the first monomer type, but any monomer type would give
      * the same answer)
      */
      VecOp::subVV(curr, system().w().rgrid(0), w0_[0]);
   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::evaluate()
   {
      system().compute();
      ++(CompressorT::mdeCounter_);
   }

   /*
   * Compute the residual for the current system state.
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::getResidual(V& resid)
   {
      // Initialize residual to -1.0
      VecOp::eqS(resid, -1.0);

      // Add c fields to get SCF residual vector
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addEqV(resid, system().c().rgrid(i));
      }
   }

   /*
   * Update the current system field coordinates.
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::update(V& newGuess)
   {
      // New field is w0_ + newGuess for the pressure field
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(wFieldTmp_[i], w0_[i], newGuess);
      }

      // Set system r-grid fields
      system().w().setRGrid(wFieldTmp_);
   }

   /*
   * Do nothing output function.
   */
   template<int D>
   void AmCompressor<D,T,V>::outputToLog()
   {}

   // Private virtual functions for vector math

   /*
   * Assign one array to another.
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::setEqual(V& a,
                                  V const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b);
   }

   /*
   * Compute and return inner product of two vectors.
   */
   template <int D, class T, class V>
   double AmCompressor<D,T,V>::dotProduct(V const & a,
                                          V const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   /*
   * Compute and return maximum absolute value element of a vector.
   */
   template <int D, class T, class V>
   double AmCompressor<D,T,V>::maxAbs(V const & a)
   {  return Reduce::maxAbs(a); }

   /*
   * Compute the vector difference a = b - c
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::subVV(V& a,
                                   V const & b,
                                   V const & c)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      VecOp::subVV(a, b, c);
   }

   /*
   * Composite a += b*c for vectors a and b, scalar c
   */
   template <int D, class T, class V>
   void AmCompressor<D,T,V>::addEqVc(V& a,
                                     V const & b,
                                     double c)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      VecOp::addEqVc(a, b, c);
   }

}
}
#endif

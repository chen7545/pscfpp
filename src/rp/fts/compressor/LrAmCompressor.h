#ifndef RP_LR_AM_COMPRESSOR_H
#define RP_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                          // base class argument
#include <pscf/iterator/AmIteratorTmpl.h>        // base class template
#include <pscf/math/IntVec.h>                    // member
#include <util/containers/DArray.h>              // member

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Anderson mixing compressor with linear-response correction step.
   *
   * Class LrAmCompressor implements an Anderson mixing algorithm in
   * which the second "correction" step is treated as a quasi-Newton
   * step, while the Jacobian is approximated by the linear response
   * of a hypothetical homogenous liquid. The residual is an r-grid
   * vector in which each element represents the deviation of the
   * sum of volume fractions from unity.
   *
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T, class V>
   class LrAmCompressor
    : public AmIteratorTmpl<typename T::Compressor, V>
   {

   public:

      /// Type of field and residual vectors.
      using VectorT = V;

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LrAmCompressor(typename T::System& system);

      /**
      * Destructor.
      */
      ~LrAmCompressor();

      /**
      * Read body of parameter file block and initialize.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. It stores the initial values of the fields
      * prior to iteration.
      *
      * \param isContinuation  true iff continuation within a sweep
      */
      void setup(bool isContinuation) override;

      /**
      * Compress to obtain partial saddle point field.
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress() override;

      /**
      * Return compressor time contributions.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const override;

      /**
      * Clear all timers and MDE solution counter.
      */
      void clearTimers() override;

   protected:

      using CompressorT = typename T::Compressor;

      // Inherited member function
      using CompressorT::system;

   private:

      /**
      * Initial values of all w fields.
      */
      DArray< typename T::RField > w0_;

      /**
      * Temporary array of w fields used in update function.
      */
      DArray< typename T::RField > wFieldTmp_;

      /**
      * Residual in real space.
      */
      typename T::RField resid_;

      /**
      * Residual in Fourier space.
      */
      typename T::RFieldDft residK_;

      /**
      * Intramolecular correlation function in Fourier space.
      */
      typename T::RField intraCorrelationK_;

      /**
      * IntraCorrelation object, used to compute intraCorrelationK_.
      */
      typename T::IntraCorrelation intra_;

      /**
      * Dimensions of wavevector mesh for real-to-complex transform.
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Number of points in k-space wavevector mesh.
      */
      int kSize_;

      /**
      * Has intraCorrelationK_ been calculated?
      */
      bool isIntraCalculated_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

      // Private AM algorithm operations

      /**
      * Compute and return the number of elements in a field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system have an initial guess for the field?
      */
      bool hasInitialGuess() override;

      /**
      * Gets the current field vector from the system.
      *
      * \param curr  current field vector
      */
      void getCurrent(VectorT& curr) override;

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate() override;

      /**
      * Compute the residual vector.
      *
      * \param resid  current residual vector value
      */
      void getResidual(VectorT& resid) override;

      /**
      * Add a correction based on the predicted residual.
      *
      * \param fieldTrial  trial field (in/out)
      * \param resTrial  predicted residual for current trial (in)
      */
      void addCorrection(VectorT& fieldTrial,
                         VectorT const & resTrial) override;

      /**
      * Update the system field with the new trial field.
      *
      * \param newGuess  trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;

      // Private virtual functions for vector math

      /**
      * Assign one field to another.
      *
      * \param a  field to be set (lhs of assignment)
      * \param b  field for it to be set to (rhs of assigment)
      */
      void setEqual(VectorT& a,
                    VectorT const & b) override;

      /**
      * Compute the inner product of two vectors
      *
      * \param a  first input vector
      * \param a  second input vector
      */
      double dotProduct(VectorT const & a,
                        VectorT const & b) override;

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(VectorT const & hist) override;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a  result vector (LHS)
      * \param b  first vector (RHS)
      * \param c  second vector (RHS)
      */
      void subVV(VectorT& a,
                 VectorT const & b,
		 VectorT const & c) override;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a  result vector (LHS)
      * \param b  input vector (RHS)
      * \param c  scalar coefficient (RHS)
      */
      void addEqVc(VectorT& a,
		   VectorT const & b,
                   double c) override;

      /// Typename alias for FFT class
      using FFTT = typename T::FFT;

      /// Typename alias for base class.
      using AmTmpl = AmIteratorTmpl< CompressorT, VectorT >;

   };

} // namespace Rp
} // namespace Pscf
#endif

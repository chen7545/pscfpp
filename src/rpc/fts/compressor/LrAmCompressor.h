#ifndef RPC_LR_AM_COMPRESSOR_H
#define RPC_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                           // base class argument
#include <pscf/iterator/AmIteratorDArray.h>       // base class template



#include <rpc/fts/compressor/IntraCorrelation.h>  // member
#include <prdc/cpu/RField.h>                      // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <pscf/math/IntVec.h>                     // member
#include <util/containers/DArray.h>               // member

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Anderson mixing compressor with linear-response correction step.
   *
   * Class LrAmCompressor implements an Anderson mixing algorithm in
   * which the second "correction" step is treated as quasi-Newton
   * step, while the Jacobian is approximated by the linear response
   * of a hypothetical homogenous liquid. The residual is an r-grid
   * vector in which each element represents the deviation of the 
   * sum of volume fractions from unity.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class LrAmCompressor
    : public AmIteratorDArray< Compressor<D> >
   {

   public:

      /// Typename for field and residual vectors.
      using VectorT = DArray<double>;

      /**
      * Constructor.
      *
      * \param system  parent System object 
      */
      LrAmCompressor(System<D>& system);

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
      * \param isContinuation true iff continuation within a sweep
      */
      void setup(bool isContinuation) override;

      /**
      * Compress to obtain partial saddle point field.
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress() override;

      /**
      * Return compressor time contribution.
      *   
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const override;

      /**
      * Clear all timers and MDE solution counter.
      */
      void clearTimers() override;

   protected:

      // Inherited protected members
      using Compressor<D>::system;

   private:

      /**
      * Initial values of all w fields.
      */
      DArray< RField<D> > w0_;

      /**
      * Temporary array of w fields used in update function.
      */
      DArray< RField<D> > wFieldTmp_;

      /**
      * Residual in real space.
      */
      RField<D> resid_;

      /**
      * Residual in Fourier space.
      */
      RFieldDft<D> residK_;

      /**
      * Intramolecular correlation function in Fourier space.
      */
      RField<D> intraCorrelationK_;

      /**
      * IntraCorrelation object, used to compute intraCorrelationK_.
      */
      IntraCorrelation<D> intra_;

      /**
      * Dimensions of wavevector mesh for real-to-complex transform.
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Number of points in k-space wavevector mesh.
      */
      int kSize_;

      /**
      * Number of times MDE has been solved for this stochastic move.
      */
      int itr_;

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
      * Add a correction based on the predicted residual.
      *
      * \param fieldTrial  trial field (in-out)
      * \param resTrial  predicted error for current trial
      */
      void addCorrection(VectorT& fieldTrial,
                         VectorT const & resTrial) override;

      /**
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess() override;

      /**
      * Gets the current field vector from the system.
      *
      * \param curr current field vector
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
      * \param resid current residual vector value
      */
      void getResidual(VectorT& resid) override;

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;

      // Indirect (grandparent) base class.
      using AmTmpl = AmIteratorTmpl<Compressor<D>, VectorT >;

   };

   // Explicit instantiation declarations
   extern template class LrAmCompressor<1>;
   extern template class LrAmCompressor<2>;
   extern template class LrAmCompressor<3>;

} // namespace Rpc
} // namespace Pscf
#endif

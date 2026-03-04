#ifndef RPC_LR_COMPRESSOR_H
#define RPC_LR_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                           // base class
#include <rpc/fts/compressor/IntraCorrelation.h>  // member
#include <prdc/cpu/RField.h>                      // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <util/containers/DArray.h>               // member
#include <util/misc/Timer.h>                      // member

namespace Pscf {
namespace Rpc {

   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Linear response compressor.
   *
   * This class implements a compressor that is an approximate Newton's
   * method in which the Jacobian is approximated by the analytically
   * calculated linear response of a homogeneous liquid with the same
   * composition as the system of interest.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class LrCompressor : public Compressor<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system System object associated with this compressor.
      */
      LrCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~LrCompressor();

      /**
      * Read all parameters and initialize.
      *
      * \param in input stream for parameter file
      */
      void readParameters(std::istream& in);

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at the
      * beginning of iteration
      */
      void setup();

      /**
      * Iterate pressure field to obtain partial saddle point.
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress();

      double subspacePercent(){return 0;};

      double correctionPercent(){return 0;};

      /**
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out) const;

      void clearTimers();

   protected:

      // Inherited protected members
      using Compressor<D>::system;
      using Compressor<D>::mdeCounter_;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::setClassName;

   private:

      // IntraCorrelation object
      IntraCorrelation<D> intra_;

      // Template w Field used in update function
      DArray< RField<D> > wFieldTmp_;

      // Residual in real space
      RField<D> resid_;

      // Residual in Fourier space
      RFieldDft<D> residK_;

      // Intramolecular correlation in Fourier space
      RField<D> intraCorrelationK_;

      // Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      // Timers for analyzing performance
      Timer timerTotal_;
      Timer timerMDE_;

      // Type of error criterion used to test convergence
      std::string errorType_;

      // Error tolerance.
      double epsilon_;

      // Current iteration counter.
      int itr_;

      // Maximum number of iterations.
      int maxItr_;

      // Total iteration counter.
      int totalItr_;

      // Verbosity level.
      int verbose_;

      // Has required memory been allocated?
      bool isAllocated_;

      // Has the IntraCorrelation been calculated?
      bool isIntraCalculated_;

      // Private member functions

      /**
      * Compute the residual vector.
      */
      void computeResidual();

      /**
      * Update system w fields.
      */
      void updateWFields();

      /**
      * Compute and return error used to test for convergence.
      *
      * \param verbose  verbosity level of output report
      * \return error  measure used to test for convergence.
      */
      double computeError(int verbose);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();

   };

   // Explicit instantiation declarations
   extern template class LrCompressor<1>;
   extern template class LrCompressor<2>;
   extern template class LrCompressor<3>;

} // namespace Rpc
} // namespace Pscf
#endif

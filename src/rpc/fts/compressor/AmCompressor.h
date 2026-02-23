#ifndef RPC_AM_COMPRESSOR_H
#define RPC_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                     // base class template argument
#include <prdc/cpu/RField.h>                // member
#include <pscf/iterator/AmIteratorDArray.h> // base class template
#include <util/containers/DArray.h>         // member



namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Anderson mixing compressor.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor
      : public AmIteratorDArray< Compressor<D> >
   {

   public:

      /// Type for state and residual vectors.
      using VectorT = DArray<double>;

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      AmCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~AmCompressor();

      /**
      * Read all parameters and initialize.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at
      * the beginning of iteration
      *
      * \param isContinuation true iff continuation within a sweep
      */
      void setup(bool isContinuation) override;

      /**
      * Compress to obtain partial saddle point w+
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress() override;

      /**
      * Return compressor times contributions.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const override;

      /**
      * Clear all timers and mde counter.
      */
      void clearTimers() override;

   protected:

      // Inherited protected member function
      using Compressor<D>::system;

   private:

      /**
      * Initial values of all fields.
      */
      DArray< RField<D> > w0_;

      /**
      * Temporary array of w fields used in update function.
      */
      DArray< RField<D> > wFieldTmp_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

      // Private virtual functions that interact with parent system

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
      * Update the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Output relevant system details to the iteration log.
      */
      void outputToLog() override;

      /// Typename alias for base class.
      using AmTmpl = AmIteratorTmpl< Compressor<D>, VectorT >;

   };

   // Explicit instantiation declarations
   extern template class AmCompressor<1>;
   extern template class AmCompressor<2>;
   extern template class AmCompressor<3>;

} // namespace Rpc
} // namespace Pscf
#endif

#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>         // memmber variable type
#include <util/containers/DArray.h>   // memmber variable type

// Forward declarations
namespace Pscf {
   namespace Correlation {
      class Mixture;
   }
   namespace Prdc {
      namespace Cpu {
         template <int D> class RField;
      }
   }
   namespace Rpc {
      template <int D> class System;
   }
}


namespace Pscf {
namespace Rpc {

   using namespace Pscf::Prdc::Cpu;

   /**
   * Intramolecular correlation analysis for LR compressors.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class IntraCorrelation
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      IntraCorrelation(System<D>& system);

      /**
      * Destructor.
      */
      ~IntraCorrelation();

      /**
      * Compute and modify intramolecular correlations.
      *
      * \param correlations  k-space grid of omega values
      */
      void computeIntraCorrelations(RField<D>& correlations);

   protected:

      /**
      * Return reference to parent system.
      */
      System<D> const & system();

   private:

      /// Pointer to a parent system object.
      System<D> const * systemPtr_;

      /// Pointer to a child Correlation::Mixture object.
      Correlation::Mixture* correlationMixturePtr_;

      /// Array of squared magnitudes for reciprocal wavevectors.
      DArray<double> Gsq_;

      /// Dimensions of Fourier grid for DFT of a real function.
      IntVec<D> kMeshDimensions_;

      /// Number of elements in the Fourier grid.
      int kSize_;

   };

   // Get the parent system by const reference.
   template <int D>
   inline System<D> const & IntraCorrelation<D>::system()
   {  return *systemPtr_; }

   // Explicit instantiation declarations
   extern template class IntraCorrelation<1>;
   extern template class IntraCorrelation<2>;
   extern template class IntraCorrelation<3>;

} // namespace Rpc
} // namespace Pscf
#endif

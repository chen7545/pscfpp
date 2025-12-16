#ifndef RPG_INTRACORRELATION_H
#define RPG_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <prdc/cuda/types.h>              // member
#include <pscf/math/IntVec.h>             // member
#include <pscf/cuda/HostDArray.h>         // member
#include <util/containers/DArray.h>       // member

// Forward references
namespace Pscf {
   namespace Correlation {
      template <typename WT> class Mixture;
   }
   namespace Prdc {
      namespace Cuda {
         template <int D> class RField; 
      }
   }   
   namespace Rpg {
      template <int D> class System;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Linear response function for response to pressure.
   *
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class IntraCorrelation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      IntraCorrelation(System<D> const & system);

      /**
      * Destructor.
      */
      ~IntraCorrelation();

      /**
      * Compute and return intramolecular correlations.
      *
      * \param correlations  k-space grid of intramolecular correlations
      */
      void computeIntraCorrelations(RField<D>& correlations);

   protected:

      /**
      * Return reference to parent system.
      */
      System<D> const & system() const;

   private:

      /// Pointer to the associated system object.
      System<D> const * systemPtr_;

      /// Pointer to a child Correlation::Mixture object.
      Correlation::Mixture<cudaReal>* correlationMixturePtr_;

      /// Array of square magnitudes for wavevectors on a k-grid.
      DArray<double> Gsq_;

      /// Host array of square magnitudes for wavevectors on a k-grid.
      HostDArray<cudaReal> correlations_;

      /// Dimensions of Fourier space grid.
      IntVec<D> kMeshDimensions_;

      /// Size of Fourier space grid. 
      int kSize_;

   };

   // Get the parent system.
   template <int D>
   inline System<D> const & IntraCorrelation<D>::system() const
   {  return *systemPtr_; }

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
#include <pscf/correlation/Mixture.h>   
namespace Pscf {
   namespace Correlation {
      extern template class Mixture<cudaReal>;
   }
   namespace Rpg {
      extern template class IntraCorrelation<1>;
      extern template class IntraCorrelation<2>;
      extern template class IntraCorrelation<3>;
   }
}
#endif

#ifndef RPG_FOURTH_ORDER_PARAMETER_H
#define RPG_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                      // base class
#include <prdc/cuda/RField.h>                     // member
#include <prdc/cuda/RFieldDft.h>                  // member
#include <pscf/math/IntVec.h>                     // member

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * FourthOrderParameter is used to detect an order-disorder transition.
   *
   * This class evaluates the sum of fourth power of the
   * Fourier mode amplitude of fluctuating fields.
   *
   * The order parameter is defined as
   * \f[
   *     \Psi_{\text{fourth}} \equiv
   *     \left[ \sum W_{-}(\bf G)^4 \right] ^{\frac{1}{4}}
   * \f]
   * where \f$W_(G)\f$ is a Fourier mode of fluctuating field.
   *
   * \see rp_FourthOrderParameter_page "Manual Page"
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class FourthOrderParameter : public AverageAnalyzer<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~FourthOrderParameter();

      /**
      * Setup before the main loop.
      */
      void setup() override;

   protected:

      /**
      * Compute and return the order parameter.
      */
      double compute() override;

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      void outputValue(int step, double value) override;

      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;

   private:

      /// Fourier transform of W_ field.
      RFieldDft<D> wK_;

      /// Prefactor for each Fourier component.
      RField<D> prefactor_;

      /// Fourth powers of Fourier magnitudes, with prefactors.
      RField<D> psi_;

      /// Dimensions of Fourier space (k-grid) mesh for real fields.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      /// Has setup been completed?
      bool isInitialized_;

      /**
      * Compute prefactor for each Fourier wavevector.
      *
      * For the real-valued function W_, each Fourier
      * coefficient G satisfies W_(G) = W_(-G). This function
      * uses Brillouin Zone (BZ) indices representation. After
      * applying fftw, if both the wavevector G and its
      * inverse -G exist in k-space, the prefactor is
      * assigned to be 1/2 for both G and -G. Otherwise,
      * it is assigned to be 1.
      */
      void computePrefactor(Array<double>& array);

      /**
      * Initialize member variable prefactor_.
      */
      void computePrefactor();

   };

   // Explicit instantiation declarations
   extern template class FourthOrderParameter<1>;
   extern template class FourthOrderParameter<2>;
   extern template class FourthOrderParameter<3>;

}
}
#endif

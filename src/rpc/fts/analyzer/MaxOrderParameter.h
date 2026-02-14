#ifndef RPC_MAX_ORDER_PARAMETER_H
#define RPC_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                      // base class
#include <prdc/cpu/RFieldDft.h>                   // member
#include <util/containers/DArray.h>               // member
#include <iostream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * MaxOrderParameter is used to detect an order-disorder transition.
   *
   * This class evalutaes maximum amplitude of the second power of the
   * Fourier mode amplitude of fluctuating fields.
   *
   * The order parameter is defined as
   * \f[
   *     \psi(k)  = \max [ |W_{-}({\bf k})|^{2} ]
   * \f]
   * where \f$ W_{-}({\bf k})\f$ is fluctuating field component with
   * wavevector \f$ {\bf k} \f$.
   *
   * \see \ref rpc_MaxOrderParameter_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class MaxOrderParameter : public AverageAnalyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      virtual ~MaxOrderParameter();

      /**
      * Setup before simulation loop.
      */
      virtual void setup();

   protected:

      /**
      * Compute and return the max order parameter.
      */
      virtual double compute();

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);

      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;

   private:

      /// W_ in Fourier mode
      RFieldDft<D> wK_;

      /// Max order parameter
      double maxOrderParameter_;

      /// q*
      IntVec<D> GminStar_;

      /// Dimensions of Fourier space (k-grid) mesh for a real field.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      /// Has setup been completed?
      bool isInitialized_;

   };

   // Explicit instantiation declarations
   extern template class MaxOrderParameter<1>;
   extern template class MaxOrderParameter<2>;
   extern template class MaxOrderParameter<3>;

}
}
#endif

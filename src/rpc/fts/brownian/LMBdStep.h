#ifndef RPC_LM_BD_STEP_H
#define RPC_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/LMBdStep.h>     // base class template
#include <rpc/system/Types.h>             // base class template argument 
#include <rpc/fts/brownian/BdStep.h>      // indirect base class
#include <prdc/cpu/RField.h>              // base class member

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class BdSimulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Leimkuhler-Mathews Brownian dynamics time stepper.
   *
   * \see Rp::LMBdStep
   * \see \ref rp_LMBdStep_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class LMBdStep : public Rp::LMBdStep<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      LMBdStep(BdSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::LMBdStep<1, Rpc::Types<1> >;
      extern template class Rp::LMBdStep<2, Rpc::Types<2> >;
      extern template class Rp::LMBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class LMBdStep<1>;
      extern template class LMBdStep<2>;
      extern template class LMBdStep<3>;
   }
}
#endif

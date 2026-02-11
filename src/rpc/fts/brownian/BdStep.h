#ifndef RPC_BD_STEP_H
#define RPC_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * This class is basically a named instantiation of the base class 
   * template Rp::BdStep, using aliases defined in Rpc::Types<D> to 
   * specialize to types used in the Rpc namespace. See documentation 
   * of the base class for details.
   *
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class BdStep : public Rp::BdStep<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator<D> object
      */
      BdStep(BdSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdStep< 1, Rpc::Types<1> >;
      extern template class BdStep< 2, Rpc::Types<2> >;
      extern template class BdStep< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class BdStep<1>;
      extern template class BdStep<2>;
      extern template class BdStep<3>;
   }
}
#endif

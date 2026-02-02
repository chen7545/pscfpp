#ifndef RPC_SIMULATOR_H
#define RPC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>    // base class template
#include <rpc/system/Types.h>              // template argument
#include <rpc/fts/simulator/SimState.h>    // member
#include <prdc/cpu/RField.h>               // member (template arg)

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Field theoretic simulator (base class).
   *
   * This class is basically a named instantiation of the base class 
   * template Rp::Simulator, using aliases defined in Rpc::Types<D> to 
   * specialize to types used in the Rpc namespace. See documentation 
   * of the base class for details.
   *
   * \see \ref rp_Simulator_page (manual page)
   *
   * \ingroup Rpc_Fts_Module
   */
   template <int D>
   class Simulator : public Pscf::Rp::Simulator<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      Simulator(typename Types<D>::System& system);

      Simulator(Simulator<D> const &) = delete;

      Simulator<D>& operator=(Simulator<D> const &) = delete;

   private:

      /// Alias for direct base class.
      using RpSimulator = Pscf::Rp::Simulator<D, Types<D> >;

   };


} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Simulator<1, Rpc::Types<1> >;
      extern template class Simulator<2, Rpc::Types<2> >;
      extern template class Simulator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Simulator<1>;
      extern template class Simulator<2>;
      extern template class Simulator<3>;
   }
}
#endif

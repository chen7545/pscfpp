#ifndef RPG_SIMULATOR_H
#define RPG_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>    // base class template
#include <rpg/system/Types.h>              // template argument
#include <rpg/fts/simulator/SimState.h>    // member
#include <prdc/cuda/RField.h>              // member (template arg)

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * Field theoretic simulator (base class).
   *
   * This class is little more than a named instantiation of the base class
   * template Rp::Simulator. It has the same public interface as the base
   * class. See the documentation of Rp::Simulator for details. 
   *
   * \ingroup Rpg_Fts_Module
   */
   template <int D>
   class Simulator : public Pscf::Rp::Simulator<D, Types<D> >
   {
   public:

      /// Alias for direct base class.
      using Base = Pscf::Rp::Simulator<D, Types<D> >;

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      Simulator(typename Types<D>::System& system);

   protected:

      /**
      * Initialize seed for vector random number generator.
      */
      void initializeVecRandom();

   };

   // Explicit instantiation declarations
   extern template class Simulator<1>;
   extern template class Simulator<2>;
   extern template class Simulator<3>;

} // namespace Rpg
} // namespace Pscf

namespace Pscf {
namespace Rp {
   // Explicit instantiation declarations for base class template
   extern template class Simulator<1, Rpg::Types<1> >;
   extern template class Simulator<2, Rpg::Types<2> >;
   extern template class Simulator<3, Rpg::Types<3> >;
}
}
#endif

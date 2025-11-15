#ifndef RPG_MIXTURE_MODIFIER_H
#define RPG_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/rl/MixtureModifier.h>  // base class template
#include "Mixture.h"                           // base class parameter

namespace Pscf {
namespace Rpg {

   /**
   * Parameter modifier for an associated Mixture.
   *
   * Class MixtureModifier<D> is derived from the template specialization 
   * Prdc::Rl::MixtureModifier< Rpg::Mixture<D> > and has the same public
   * interface as this base class. See documentation of this base class
   * template for details.
   *
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class MixtureModifier : public Rl::MixtureModifier< Mixture<D> >
   {

   public:
 
      /// Direct (parent) base class.
      using Base = typename Prdc::Rl::MixtureModifier< Mixture<D> >;

      // Inherited public member functions

      using Base::associate;
      using Base::setKuhn;
      using Base::setPhiPolymer;
      using Base::setMuPolymer;
      using Base::setBlockLength;
      using Base::setPhiSolvent;
      using Base::setMuSolvent;
      using Base::setSolventSize;
      using Base::setVMonomer;

   };

   // Explicit instantiation declarations
   extern template class MixtureModifier<1>;
   extern template class MixtureModifier<2>;
   extern template class MixtureModifier<3>;

} // namespace Rpg

namespace Prdc {
   // Explicit instantiation declarations for base class 
   extern template class Rl::MixtureModifier< Rpg::Mixture<1> >;
   extern template class Rl::MixtureModifier< Rpg::Mixture<2> >;
   extern template class Rl::MixtureModifier< Rpg::Mixture<3> >;
} // namespace Prdc

} // namespace Pscf
#endif

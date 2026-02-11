#ifndef RPG_ITERATOR_TPP
#define RPG_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/environment/Environment.h>

#include <rp/scft/iterator/Iterator.tpp>

namespace Pscf {
namespace Rpg {

   /*
   * Default constructor.
   */
   template <int D>
   Iterator<D>::Iterator()
    : Rp::Iterator<D, System<D> >()
   {}

   /*
   * Constructor.
   */
   template <int D>
   Iterator<D>::Iterator(System<D>& system)
    : Rp::Iterator<D, System<D> >(system)
   {}

} // namespace Rpg
} // namespace Pscf
#endif

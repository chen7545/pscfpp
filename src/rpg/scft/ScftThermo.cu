/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"            // class header
#include <rp/scft/ScftThermo.tpp>  // base class template implementation

#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>
#include <pscf/cuda/Reduce.h>
#include <prdc/cuda/RField.h>
#include <pscf/interaction/Interaction.h>

namespace Pscf {

   namespace Rp {
      // Explicit instantiation definitions for base class
      template class ScftThermo<1, Rpg::System<1> >;
      template class ScftThermo<2, Rpg::System<2> >;
      template class ScftThermo<3, Rpg::System<3> >;
   }

   namespace Rpg {

      /*
      * Constructor
      */
      template <int D>
      ScftThermo<D>::ScftThermo(System<D> const & system)
       : Base(system)
      {};

      /*
      * Inner product of r-grid fields.
      */
      template <int D>
      double ScftThermo<D>::innerProduct(RFieldT const & A, 
                                         RFieldT const & B) const
      {  return Cuda::Reduce::innerProduct(A, B); };

      // Explicit instantiation
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;

   }
}

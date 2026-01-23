/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"            // class header
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/ScftThermo.tpp>  // base class implementation

namespace Pscf {

   namespace Rp {
      // Explicit instantiation of base class template
      template class ScftThermo<1, Rpc::System<1> >;
      template class ScftThermo<2, Rpc::System<2> >;
      template class ScftThermo<3, Rpc::System<3> >;
   }

   namespace Rpc {

      /*
      * Constructor
      */
      template <int D>
      ScftThermo<D>::ScftThermo(System<D> const & system)
       : Base(system)
      {};

      #if 0
      /*
      * Inner product of r-grid fields.
      */
      template <int D>
      double ScftThermo<D>::innerProduct(RFieldT const & A,
                                         RFieldT const & B) const
      {
         const int meshSize = Base::domain().mesh().size();
         UTIL_CHECK(meshSize == A.capacity())
         UTIL_CHECK(meshSize == B.capacity())
         double sum = 0.0;
         for (int k = 0; k < meshSize; ++k) {
            sum += A[k]*B[k];
         }
         return sum;
      };
      #endif

      // Explicit instantiation
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;

   }
}

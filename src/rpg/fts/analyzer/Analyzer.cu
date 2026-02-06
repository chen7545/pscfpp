/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/Analyzer.tpp>

// Subclass constructor definition
namespace Pscf {
   namespace Rpg {

      /*
      * Constructor.
      */
      template <int D>
      Analyzer<D>::Analyzer(Simulator<D>& simulator, System<D>& system)
       : Rp::Analyzer<D, Simulator<D>, System<D> >(simulator, system)
      {}

   }
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class Analyzer<1, Rpg::Simulator<1>, Rpg::System<1> >;
      template class Analyzer<2, Rpg::Simulator<2>, Rpg::System<2> >;
      template class Analyzer<3, Rpg::Simulator<3>, Rpg::System<3> >;
   } 
   namespace Rpg {
      template class Analyzer<1>;
      template class Analyzer<2>;
      template class Analyzer<3>;
   } 
} 

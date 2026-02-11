#ifndef RPC_HAMILTONIAN_ANALYZER_H
#define RPC_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"                  // indirect base class
#include <rp/fts/analyzer/HamiltonianAnalyzer.h>  // base class template
#include <rpc/system/Types.h>                     // template argument

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Compute averages and output block averages of Hamiltonian components.
   *
   * Instantiations of this template are basically named instantiations 
   * of the base class template Rp::HamiltonianAnalyzer, using type name
   * aliases defined by the Rpc::Types<D> class. See the documentation 
   * for this base class template for details. 
   *
   * \see \ref rp_HamiltonianAnalyzer_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class HamiltonianAnalyzer 
     : public Rp::HamiltonianAnalyzer< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      HamiltonianAnalyzer(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class HamiltonianAnalyzer< 1, Rpc::Types<1> >;
      extern template class HamiltonianAnalyzer< 2, Rpc::Types<2> >;
      extern template class HamiltonianAnalyzer< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class HamiltonianAnalyzer<1>;
      extern template class HamiltonianAnalyzer<2>;
      extern template class HamiltonianAnalyzer<3>;
   }
}
#endif

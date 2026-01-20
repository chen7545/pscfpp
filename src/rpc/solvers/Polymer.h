#ifndef RPC_POLYMER_H
#define RPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <rpc/system/Types.h>      // base class template

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Descriptor and solver for one polymer species.
   *
   * This class is simply a named instantiation of the base class 
   * template Rp::Polymer. Please see documentation of the base class
   * for details.
   *
   * \ref user_param_polymer_sec "Manual Page"
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Polymer : public Rp::Polymer<D, Types<D>>
   {

   public:


      /// Base class, partial template specialization.
      using Base = PolymerTmpl< Block<D>, Propagator<D> >;

      #if 0
      // Inherited public member functions
      using Base::edge;
      using Base::block;
      using Base::propagator;
      using PolymerSpecies<double>::vertex;
      using PolymerSpecies<double>::propagatorId;
      using PolymerSpecies<double>::path;
      using PolymerSpecies<double>::nBlock;
      using PolymerSpecies<double>::nVertex;
      using PolymerSpecies<double>::nPropagator;
      using PolymerSpecies<double>::length;
      using PolymerSpecies<double>::nBead;
      using PolymerSpecies<double>::type;
      using Species<double>::phi;
      using Species<double>::mu;
      using Species<double>::q;
      using Species<double>::ensemble;
      using Species<double>::setPhi;
      using Species<double>::setMu;
      #endif

   };

   // Explicit instantiation declarations
   extern template class Polymer<1>;
   extern template class Polymer<2>;
   extern template class Polymer<3>;

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations for base classes
namespace Pscf {
   namespace Rp {
      extern template class Polymer<1, Rpc::Types<1> >;
      extern template class Polymer<2, Rpc::Types<2> >;
      extern template class Polymer<3, Rpc::Types<3> >;
   } 
}

#endif

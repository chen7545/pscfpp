#ifndef RPC_PROPAGATOR_H
#define RPC_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Propagator.h>  // base class template
#include <rpc/system/Types.h>       // base class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * MDE solver for one direction of one block.
   *
   * Most of the functionality of this class is inherited from the base
   * class Rp::Propagator<D, Rpc::Types<D> >, in which Types<D> is a
   * class that contains a set of typename aliases for use in the Rpc
   * program level namespace. See the documentation of the Rp::Propagator 
   * class template for details. 
   *
   * Virtual allocate() and reallocate() functions are defined here and 
   * in the analogous class Rpg::Propagator<D> because the CPU ad GPU
   * versions of this class use different strategies for memory allocation.
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Propagator : public Rp::Propagator< D, Types<D> >
   {

   public:

      // Member functions

      /**
      * Constructor.
      */
      Propagator();

      /**
      * Destructor.
      */
      ~Propagator();

      /**
      * Allocate memory used by this propagator.
      *
      * The parameter ns is the number of values of s at which q(r,s) is
      * calculated, including the end values at the terminating vertices.
      * See docs for the function ns(), which returns this value.
      *
      * The address of the associated Mesh<D> object is retained.
      *
      * An Exception is thrown if the propagator is already allocated.
      *
      * \param ns  number of slices (including end points at vertices)
      * \param mesh  spatial discretization mesh
      */
      void allocate(int ns, const Mesh<D>& mesh) override;

      /**
      * Reallocate memory used by this propagator.
      *
      * This function is used when the value of ns is changed, which can
      * occur during some parameter sweeps. See docs for allocate and ns.
      *
      * An Exception is thrown if the propagator has not been previously
      * allocated, or if it is allocated but the value of ns is unchanged.
      *
      * \param ns  number of slices (including end points)
      */
      void reallocate(int ns) override;

   protected:

      /// Direct (parent) base class.
      using RpPropagatorT = Rp::Propagator<D, Types<D> >;

      /// Inherited typename alias for indirect (grandparent) base class.
      using typename RpPropagatorT::PropagatorTmplT;

      /// Inherited non-dependent class members.
      using RpPropagatorT::qFields_;
      using RpPropagatorT::ns_;
      using RpPropagatorT::isAllocated_;
      using RpPropagatorT::computeHead;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Propagator<1, Rpc::Types<1> >;
      extern template class Propagator<2, Rpc::Types<2> >;
      extern template class Propagator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Propagator<1>;
      extern template class Propagator<2>;
      extern template class Propagator<3>;
   }
}
#endif

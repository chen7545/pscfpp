#ifndef RP_PROPAGATOR_H
#define RP_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <util/containers/DArray.h>      // member


// Forward declaration
namespace Pscf {
   template <int D, class T> class Mesh;
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * MDE solver for one direction of one block.
   *
   * A fully initialized Propagator<D> has an associations with a Block
   * object that owns this propagator and its partner, and with a partner
   * Propagator<D> that solves the MDE within the same block in the
   * opposite direction. It also has an association with a Mesh<D> that
   * describes a spatial grid, and associations with zero or more source
   * Propagator<D> objects that are used to compute an initial condition
   * for this propagator at the head vertex.
   *
   * The associated Block stores information required to numerically
   * solve the modified diffusion equation (MDE), including quantities
   * that depend upon the w-field associated with this block, the unit
   * cell parameters and (in the thread model) the contour step size.
   * These quantities are set and stored by the block because their values
   * are the same for the two propagators owned by each block, but may be
   * different for different blocks.  The algorithm used by a Propagator
   * to solve the MDE repeatedly calls step functions provided by the
   * parent Block.
   *
   * \ingroup Rp_Solver_Module
   */
   template <int D, class T>
   class Propagator : public PropagatorTmpl< T::Propagator >
   {

   public:

      // Public typename aliases

      // Direct (parent) base class.
      using PropagatorTmplT = PropagatorTmpl< T::Propagator >;

      // Types specific to program-level namespace
      using FieldT = T::RField;
      using BlockT = T::Block;

      // Member functions

      /**
      * Constructor.
      */
      Propagator();

      /**
      * Destructor.
      */
      virtual ~Propagator();

      /**
      * Associate this propagator with a block.
      *
      * \param block associated Block object
      */
      void setBlock(BlockT& block);

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
      virtual void allocate(int ns, const Mesh<D>& mesh);

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
      virtual void reallocate(int ns);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial q-field at the head of this
      * block, and then solves the modified diffusion equation (MDE) to
      * propagate the solution from the head to the tail. Algorithms for
      * the thread or bead model may be used, depending on the value of
      * PolymerModel::model().
      */
      void solve();

      /**
      * Solve the MDE for a specified initial condition.
      *
      * This function solves the modified diffusion equation (MDE) for
      * this block with a specified initial condition, which is given by
      * the function parameter "head".
      *
      * \param head  initial condition of q-field at head of block
      */
      void solve(FieldT const & head);

      /**
      * Compute and return partition function for the polymer molecule.
      *
      * This function computes the partition function Q for the molecule
      * as a spatial average of the pointwise product of the initial/head
      * slice for this propagator and the final/tail slice of its partner.
      *
      * \param Q  output value, spatial average of q*q^{+} at head
      */
      void computeQ(double & Q) const;

      /**
      * Return q-field at a specified step.
      *
      * \param i  step index, 0 <= i < ns
      */
      const FieldT& q(int i) const;

      /**
      * Return q-field at the initial (head) vertex.
      */
      const FieldT& head() const;

      /**
      * Return q-field at the terminal (tail) vertex.
      *
      * This function throws an Exception if invoked while the bead model
      * is in use (i.e., if PolymerModel::isThread() == false) and the tail
      * for this propagator is a chain end (i.e., if isTailEnd() == true).
      * In this case, the tail slice is not needed, and so is not computed.
      */
      const FieldT& tail() const;

      /**
      * Get the associated Block object by const reference.
      */
      BlockT const & block() const;

      /**
      * Get the number of values of s (or slices), including head and tail.
      *
      * The value of ns is the number of values of s at which q(r,s) is
      * calculated, including the end values at the terminating vertices
      * (the head and tail).  In the bead model, this is two more than the
      * number of beads in the block. In the thread model, this is one
      * more than the number length/ds of contour length steps.
      */
      int ns() const;

      /**
      * Has memory been allocated for this propagator?
      */
      bool isAllocated() const;

      // Inherited public members with non-dependent names

      using PropagatorTmplT::nSource;
      using PropagatorTmplT::source;
      using PropagatorTmplT::partner;
      using PropagatorTmplT::setIsSolved;
      using PropagatorTmplT::isSolved;
      using PropagatorTmplT::hasPartner;
      using PropagatorTmplT::isHeadEnd;
      using PropagatorTmplT::isTailEnd;

   protected:

      /// Array of propagator slices at different contour variable values.
      DArray<FieldT> qFields_;

      /// Number of slices, including head and tail slices.
      int ns_;

      /// Is this propagator allocated?
      bool isAllocated_;

      /**
      * Compute initial q-field at the head vertex.
      *
      * In either the thread or bead model, the head slice of each
      * propagator is the product of tail slices for incoming propagators
      * from other bonds that terminate at the head vertex, or 1 for a
      * chain end.
      */
      void computeHead();

   private:

      /// Pointer to the associated Block.
      BlockT* blockPtr_;

      /// Pointer to the associated Mesh.
      Mesh<D> const * meshPtr_;

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   template <int D, class T> inline
   typename Propagator<D,T>::FieldT const& Propagator<D,T>::head() const
   {
      UTIL_CHECK(isSolved());
      return qFields_[0];
   }

   /*
   * Return q-field at end of block.
   */
   template <int D, class T> inline
   typename Propagator<D,T>::FieldT const& Propagator<D,T>::tail() const
   {
      UTIL_CHECK(isSolved());
      UTIL_CHECK(PolymerModel::isThread() || !isTailEnd());
      return qFields_[ns_-1];
   }

   /*
   * Return q-field at specified step.
   */
   template <int D, class T> inline
   typename Propagator<D,T>::FieldT const& Propagator<D,T>::q(int i) const
   {
      UTIL_CHECK(isSolved());
      return qFields_[i];
   }

   /*
   * Get the associated Block object by const reference.
   */
   template <int D, class T> inline
   typename T::Block const & Propagator<D,T>::block() const
   {
      UTIL_ASSERT(blockPtr_);
      return *blockPtr_;
   }

   /*
   * Get the number of counter grid points.
   */
   template <int D, class T> inline
   int Propagator<D,T>::ns() const
   {  return ns_; }

   /*
   * Has memory been allocated for this propagator?
   */
   template <int D, class T> inline
   bool Propagator<D,T>::isAllocated() const
   {  return isAllocated_; }

   /*
   * Associate this propagator with a unique block.
   */
   template <int D, class T> inline
   void Propagator<D,T>::setBlock(BlockT& block)
   {  blockPtr_ = &block; }

} // namespace Rp
} // namespace Pscf
#endif

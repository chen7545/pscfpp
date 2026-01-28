#ifndef RP_PROPAGATOR_TPP
#define RP_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   Propagator<D,T>::Propagator()
    : ns_(0),
      isAllocated_(false),
      blockPtr_(nullptr),
      meshPtr_(nullptr)
   {}

   /*
   * Destructor.
   */
   template <int D, class T>
   Propagator<D,T>::~Propagator()
   {}

   /*
   * Allocate memory used by this propagator.
   */
   template <int D, class T>
   void Propagator<D,T>::allocate(int ns, const Mesh<D>& mesh)
   {
      UTIL_CHECK(!isAllocated_);
      ns_ = ns;
      meshPtr_ = &mesh;
   }

   /*
   * Reallocate memory used by this propagator using new ns value.
   */
   template <int D, class T>
   void Propagator<D,T>::reallocate(int ns)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(ns_ != ns);
      ns_ = ns;
   }

   /*
   * Compute initial head q-field for the thread model.
   */
   template <int D, class T>
   void Propagator<D,T>::computeHead()
   {
      UTIL_CHECK(meshPtr_);

      // Reference to head slice of this propagator
      typename T::RField& qh = qFields_[0];

      // Initialize head slice qh to 1.0 at all grid points
      VecOp::eqS(qh, 1.0);

      // Pointwise multiply tail q-fields of all sources
      if (!isHeadEnd()) {
         const int ns = nSource();
         UTIL_CHECK(ns > 0);
         for (int is = 0; is < ns; ++is) {
            if (!source(is).isSolved()) {
               UTIL_THROW("Source not solved in computeHead");
            }
            VecOp::mulEqV(qh, source(is).tail());
         }
      }

   }

   /*
   * Solve the modified diffusion equation for this block.
   */
   template <int D, class T>
   void Propagator<D,T>::solve()
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(isAllocated());

      // Initialize head as pointwise product of source propagators
      computeHead();

      if (PolymerModel::isThread()) {

         // MDE step loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else 
      if (PolymerModel::isBead()) {

         // Half-bond and bead weight for first bead
         if (isHeadEnd()) {
            VecOp::eqV(qFields_[1], qFields_[0]);
         } else {
            block().stepHalfBondBead(qFields_[0], qFields_[1]);
         }
         block().stepFieldBead(qFields_[1]);

         // MDE step loop for bead model (stop before tail vertex)
         int iStep;
         for (iStep = 1; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Half-bond for tail slice
         if (isTailEnd()) {
            VecOp::eqV(qFields_[ns_-1], qFields_[ns_-2]);
         } else {
            block().stepHalfBondBead(qFields_[ns_-2], qFields_[ns_-1]);
         }

      } else {
         // This should be impossible
         UTIL_THROW("Unexpected PolymerModel type");
      }

      PropagatorTmplT::setIsSolved(true);
   }

   /*
   * Solve the MDE with a specified initial condition at the head.
   */
   template <int D, class T>
   void Propagator<D,T>::solve(typename T::RField const & head)
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(head.capacity() == mesh().size());

      // Initialize initial (head) slice
      VecOp::eqV(qFields_[0], head);

      if (PolymerModel::isThread()) {

         // MDE integration loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else 
      if (PolymerModel::isBead()) {

         // Half-bond and bead weight for first bead
         if (isHeadEnd()) {
            VecOp::eqV(qFields_[1], qFields_[0]);
         } else {
            block().stepHalfBondBead(qFields_[0], qFields_[1]);
         }
         block().stepFieldBead(qFields_[1]);

         // MDE step loop for bead model (stop before tail vertex)
         int iStep;
         for (iStep = 1; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Half-bond for tail
         if (isTailEnd()) {
            VecOp::eqV(qFields_[ns_-1], qFields_[ns_-2]);
         } else {
            block().stepHalfBondBead(qFields_[ns_-2], qFields_[ns_-1]);
         }

      } else {
         // This should be impossible
         UTIL_THROW("Unexpected PolymerModel type");
      }

      PropagatorTmplT::setIsSolved(true);
   }

   /*
   * Compute spatial average of product of head and tail of partner.
   */
   template <int D, class T>
   void Propagator<D,T>::computeQ(double & Q) const
   {
      // Preconditions
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(isAllocated_);
      if (!isSolved()) {
         UTIL_THROW("Propagator is not solved.");
      }
      if (!hasPartner()) {
         UTIL_THROW("Propagator has no partner set.");
      }
      if (!partner().isSolved()) {
         UTIL_THROW("Partner propagator is not solved");
      }
      UTIL_CHECK(isHeadEnd() == partner().isTailEnd());

      Q = 0.0;
      if (PolymerModel::isBead() && isHeadEnd()) {
         // Compute average of q for last bead of partner
         Q = Reduce::sum(partner().q(ns_-2)); 
      } else {
         // Compute average product of head slice and partner tail slice
         Q = Reduce::innerProduct(head(), partner().tail()); 
      }
      Q /= double(mesh().size());
   }

} // namespace Rp
} // namespace Pscf
#endif

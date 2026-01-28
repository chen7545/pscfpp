#ifndef RP_PROPAGATOR_TPP
#define RP_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   Propagator<D,T>::Propagator()
    : blockPtr_(nullptr),
      meshPtr_(nullptr),
      ns_(0),
      isAllocated_(false),
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
      FieldT& qh = qFields_[0];

      // Initialize head slice qh to 1.0 at all grid points
      VecOp::eqS(qh, 1.0);

      // Pointwise multiply tail q-fields of all sources
      if (!isHeadEnd()) {
         UTIL_CHECK(nSource() > 0);
         for (int is = 0; is < nSource(); ++is) {
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

      setIsSolved(true);
   }

   /*
   * Solve the MDE with a specified initial condition at the head.
   */
   template <int D, class T>
   void Propagator<D,T>::solve(FieldT const & head)
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(meshPtr_);
      int nx = meshPtr_->size();
      UTIL_CHECK(head.capacity() == nx);

      // Initialize initial (head) field
      FieldT& qh = qFields_[0];
      for (int i = 0; i < nx; ++i) {
         qh[i] = head[i];
      }

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

      setIsSolved(true);
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
      //int nx = meshPtr_->size();

      Q = 0.0;
      if (PolymerModel::isBead() && isHeadEnd()) {
         // Compute average of q for last bead of partner
         FieldT const& qt = partner().q(ns_-2);
         Q = Reduce::sum(qt); 
      } else {
         // Compute average product of head slice and partner tail slice
         FieldT const& qh = head();
         FieldT const& qt = partner().tail();
         Q = Reduce::innerProduct(qh, qt); 
      }
      Q /= double(meshPtr_->size());
   }

} // namespace Rp
} // namespace Pscf
#endif

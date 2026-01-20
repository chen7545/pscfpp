#ifndef RPG_SOLVENT_TPP
#define RPG_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <prdc/cuda/resources.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rpg { 

   /*
   * Constructor.
   */
   template <int D>
   Solvent<D>::Solvent()
    : meshPtr_(nullptr)
   {  ParamComposite::setClassName("Solvent"); }

   /*
   * Destructor.
   */
   template <int D>
   Solvent<D>::~Solvent()
   {}

   /*
   * Create an association with a mesh.
   */
   template <int D> 
   void Solvent<D>::associate(Mesh<D> const & mesh)
   {
      UTIL_CHECK(!meshPtr_)
      UTIL_CHECK(mesh.size() > 1);
      meshPtr_ = &mesh;
   }

   /*
   * Allocate memory for the concentration field (cField).
   */
   template <int D>
   void Solvent<D>::allocate()
   {
      UTIL_CHECK(meshPtr_);
      cField_.allocate(meshPtr_->dimensions());
   }

   /*
   * Compute concentration, q, and phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(RField<D> const & wField, double phiTot)
   {
      // Local constants
      const int nx = meshPtr_->size();
      const double s = SolventSpecies<cudaReal>::size();

      // Evaluate unnormalized integral and Q
      double Q = 0.0;
      VecOp::expVc(cField_, wField, -1.0*size());
      Q = Reduce::sum(cField_);
      Q = Q / double(nx);     // spatial average
      Q /= phiTot;            // correct for partial occupation

      // Note: phiTot = 1.0 except in the case of a mask that confines
      // material to a fraction of the unit cell. 

      // Set q and compute mu or phi
      Species::setQ(Q);

      // Normalize concentration 
      double prefactor = phi()/Q;
      VecOp::mulEqS(cField_, prefactor);
   }

}
}
#endif

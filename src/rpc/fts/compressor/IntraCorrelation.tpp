#ifndef RPC_INTRACORRELATION_TPP
#define RPC_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/correlation/Mixture.h>

#include <util/global.h>

namespace Pscf {
namespace Rpc{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   IntraCorrelation<D>::IntraCorrelation(System<D>& system)
    : systemPtr_(&system),
      correlationMixturePtr_(nullptr),
      kSize_(-1)
   {
      correlationMixturePtr_ 
          = new Pscf::Correlation::Mixture<double>(system.mixture());
   }

   /*
   * Destructor.
   */
   template <int D>
   IntraCorrelation<D>::~IntraCorrelation()
   {
      delete correlationMixturePtr_;
   }

   /*
   * Compute k-space array of intramolecular correlation functions.
   */
   template<int D>
   void
   IntraCorrelation<D>::computeIntraCorrelations(RField<D>& correlations)
   {
      // Local copies of domain properties
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Compute Fourier space kMeshDimensions_ and kSize_
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);
      UTIL_CHECK(correlations.capacity() == kSize_);

      // Check allocation of Gsq_ (k-space array of square wavenumbers)
      if (!Gsq_.isAllocated()) {
         Gsq_.allocate(kSize_);
      }
      UTIL_CHECK(Gsq_.capacity() == kSize_);

      // Compute Gsq_
      IntVec<D> G, Gmin;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, dimensions, unitCell);
         Gsq_[iter.rank()] = unitCell.ksq(Gmin);
      }

      // Compute total intramolecular correlation function
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);

   }

}
}
#endif

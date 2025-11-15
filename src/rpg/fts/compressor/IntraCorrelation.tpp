#ifndef RPG_INTRACORRELATION_TPP
#define RPG_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/types.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/SolventSpecies.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/EdgeIterator.h>
#include <pscf/correlation/Debye.h>
#include <pscf/correlation/Mixture.h>

#include <pscf/cuda/HostDArray.h>

#include <util/global.h>

namespace Pscf {
namespace Rpg{

   using namespace Util;
   using namespace Prdc; 
   using namespace Prdc::Cuda; 

   /*
   * Constructor.
   */
   template <int D>
   IntraCorrelation<D>::IntraCorrelation(System<D> const & system)
    : systemPtr_(&system),
      correlationMixturePtr_(nullptr),
      kSize_(-1)
   {
      setClassName("IntraCorrelation");
      correlationMixturePtr_ = new Correlation::Mixture<cudaReal>(system.mixture());
   }

   /*
   * Destructor.
   */
   template <int D>
   IntraCorrelation<D>::~IntraCorrelation()
   {  delete correlationMixturePtr_; }

   template<int D>
   void
   IntraCorrelation<D>::computeIntraCorrelations(RField<D>& correlations)
   {
      // Local copies of system properties
      Mixture<D> const & mixture = system().mixture();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Compute k-space mesh dimensions kMeshDimensions_ & size kSize_
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);
      UTIL_CHECK(correlations.capacity() == kSize_);

      // Check allocation of Gsq_ (k-space array of square wavenumber)
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

      // Check allocation of HostDArray<double> correlations_
      if (!correlations_.isAllocated()) {
         correlations_.allocate(kSize_);
      }
      UTIL_CHECK(correlations_.capacity() == kSize_);

      // Compute array of correlation function values on CPU
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations_);

      #if 0
      // Initialize correlations_ to zero
      for (int i = 0; i < correlations_.capacity(); ++i){
         correlations_[i] = 0.0;
      }

      const double vMonomer = mixture.vMonomer();
      const int nPolymer = mixture.nPolymer();
      const int nSolvent = mixture.nSolvent();

      double phi, cPolymer, polymerLength, rsqAB;
      double length, lengthA, lengthB, lengthC, ksq;
      double kuhn, kuhnA, kuhnB, kuhnC, eA, eB;
      int monomerId, monomerIdA, monomerIdB, monomerIdC;
      int rank;

      // Loop over polymer species
      for (int i = 0; i < nPolymer; i++){

         // Local copies of polymer properties
         PolymerSpecies const & polymer = mixture.polymerSpecies(i);
         const int nBlock = polymer.nBlock();
         phi = polymer.phi();

         // Compute polymerLength (sum of lengths of all blocks)
         polymerLength = 0.0;
         for (int j = 0; j < nBlock; j++) {
            if (PolymerModel::isThread()) {
               length = polymer.edge(j).length();
            } else {
               length = (double) polymer.edge(j).nBead();
            }
            polymerLength += length;
         }

         // Compute cPolymer (polymer number concentration)
         cPolymer = phi/(polymerLength*vMonomer);

         // Compute diagonal (j = k) contributions
         for (int j = 0; j < nBlock; j++) {

            monomerId = polymer.edge(j).monomerId();
            kuhn = mixture.monomer(monomerId).kuhn();

            // Loop over ksq to increment correlations_
            if (PolymerModel::isThread()) {
               length = polymer.edge(j).length();
               for (iter.begin(); !iter.atEnd(); ++iter) {
                  rank = iter.rank();
                  ksq = Gsq_[rank];
                  correlations_[rank] +=
                            cPolymer * Correlation::dt(ksq, length, kuhn);
               }
            } else {
               length = (double) polymer.edge(j).nBead();
               for (iter.begin(); !iter.atEnd(); ++iter) {
                  rank = iter.rank();
                  ksq = Gsq_[rank];
                  correlations_[rank] +=
                            cPolymer * Correlation::db(ksq, length, kuhn);
               }
            }

         }

         // Compute off-diagonal contributions
         if (nBlock > 1) {
            EdgeIterator edgeItr(polymer);

            // Outer loop over blocks
            for (int ia = 1; ia < nBlock; ++ia) {

               // Block A properties
               Edge const & edgeA = polymer.edge(ia);
               if (PolymerModel::isThread()) {
                  lengthA = edgeA.length();
               } else {
                  lengthA = (double) edgeA.nBead();
               }
               monomerIdA = edgeA.monomerId();
               kuhnA = mixture.monomer(monomerIdA).kuhn();

               // Inner loop over blocks
               for (int ib = 0; ib < ia; ++ib)  {

                  // Block B properties
                  Edge const & edgeB = polymer.edge(ib);
                  if (PolymerModel::isThread()) {
                     lengthB = edgeB.length();
                  } else {
                     lengthB = (double) edgeB.nBead();
                  }
                  monomerIdB = edgeB.monomerId();
                  kuhnB = mixture.monomer(monomerIdB).kuhn();

                  // Initialize rsqAB
                  if (PolymerModel::isThread()) {
                     rsqAB = 0.0;
                  } else {
                     rsqAB = 0.5*(kuhnA*kuhnA + kuhnB*kuhnB);
                  }

                  // Iterate over intermediate blocks, if any
                  int edgeId;
                  edgeItr.begin(ia, ib);
                  while (edgeItr.notEnd()) {
                     edgeId = edgeItr.currentEdgeId();
                     Edge const & edgeC = polymer.edge(edgeId);
                     monomerIdC = edgeC.monomerId();
                     kuhnC = mixture.monomer(monomerIdC).kuhn();
                     if (edgeId != ia && edgeId != ib) {
                        if (PolymerModel::isThread()) {
                           lengthC = edgeC.length();
                        } else {
                           lengthC = double(edgeC.nBead());
                        }
                        rsqAB += lengthC * kuhnC * kuhnC;
                     }
                     ++edgeItr;
                  }

                  // Loop over ksq to increment intra-block correlations
                  double x;
                  double prefactor = 2.0*cPolymer;
                  if (PolymerModel::isThread()) {
                     for (iter.begin(); !iter.atEnd(); ++iter) {
                        rank = iter.rank();
                        ksq = Gsq_[rank];
                        x = std::exp( -rsqAB * ksq / 6.0);
                        eA = Correlation::et(ksq, lengthA, kuhnA);
                        eB = Correlation::et(ksq, lengthB, kuhnB);
                        correlations_[rank] += prefactor * x * eA * eB;
                     }
                  } else {
                     for (iter.begin(); !iter.atEnd(); ++iter) {
                        rank = iter.rank();
                        ksq = Gsq_[rank];
                        x = std::exp( -rsqAB * ksq / 6.0);
                        eA = Correlation::eb(ksq, lengthA, kuhnA);
                        eB = Correlation::eb(ksq, lengthB, kuhnB);
                        correlations_[rank] += prefactor * x * eA * eB;
                     }
                  }

               } // loop: block ib
            } // loop: block ia
         } // if (nBlock > 1) - off diagonal elements

      } // loop over polymer species

      // Loop over solvent species (if any)
      if (nSolvent > 0) {
         double size, dcorr;
         for (int i = 0; i < nSolvent; i++){
            SolventSpecies const & solvent = mixture.solventSpecies(i);
            phi = solvent.phi();
            size = solvent.size();
            dcorr = phi*size/vMonomer;
            for (iter.begin(); !iter.atEnd(); ++iter) {
               correlations_[iter.rank()] += dcorr;
            }
         }
      }
      #endif

      // Copy host array "correlations_" to device array "correlations"
      correlations = correlations_;

   }

}
}
#endif

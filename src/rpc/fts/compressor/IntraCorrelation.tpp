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
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/SolventSpecies.h>
#include <pscf/chem/Edge.h>
#include <pscf/correlation/Debye.h>
#include <pscf/correlation/Mixture.h>
#include <pscf/chem/EdgeIterator.h>

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
      correlationMixturePtr_ = new Correlation::Mixture(system.mixture());
   }

   /*
   * Destructor.
   */
   template <int D>
   IntraCorrelation<D>::~IntraCorrelation()
   {}

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

      // Compute total Omega
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);

      #if 0
      // Initialize correlations to zero
      for (int i = 0; i < correlations.capacity(); ++i){
         correlations[i] = 0.0;
      }

      Mixture<D> const & mixture = system().mixture();
      const int nPolymer = mixture.nPolymer();
      const int nSolvent = mixture.nSolvent();
      const double vMonomer = mixture.vMonomer();

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

            // Loop over ksq to increment correlations
            if (PolymerModel::isThread()) {
               length = polymer.edge(j).length();
               for (iter.begin(); !iter.atEnd(); ++iter) {
                  rank = iter.rank();
                  ksq = Gsq_[rank];
                  correlations[rank] +=
                             cPolymer * Correlation::dt(ksq, length, kuhn);
               }
            } else {
               length = (double) polymer.edge(j).nBead();
               for (iter.begin(); !iter.atEnd(); ++iter) {
                  rank = iter.rank();
                  ksq = Gsq_[rank];
                  correlations[rank] +=
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
                        correlations[rank] += prefactor * x * eA * eB;
                     }
                  } else {
                     for (iter.begin(); !iter.atEnd(); ++iter) {
                        rank = iter.rank();
                        ksq = Gsq_[rank];
                        x = std::exp( -rsqAB * ksq / 6.0);
                        eA = Correlation::eb(ksq, lengthA, kuhnA);
                        eB = Correlation::eb(ksq, lengthB, kuhnB);
                        correlations[rank] += prefactor * x * eA * eB;
                     }
                  }

               } // loop: block ib
            } // loop: block ia
         } // if (nBlock > 1)

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
               correlations[iter.rank()] += dcorr;
            }
         }
      }
      #endif

      #if 0
      //DArray<double> corrTest;
      //corrTest.allocate(kSize_);
      //correlationMixturePtr_->computeOmegaTotal(Gsq_, corrTest);

      for (int i = 0; i < kSize_; ++i) {
         UTIL_CHECK(std::abs(corrTest[i] - correlations[i]) < 1.0E-8);
      }
      #endif

   }

}
}
#endif

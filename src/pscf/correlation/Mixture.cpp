#ifndef PSCF_OMEGA_MIXTURE_CPP
#define PSCF_OMEGA_MIXTURE_CPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include "Polymer.h"

#include <pscf/correlation/Debye.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/SolventSpecies.h>
#include <pscf/chem/MixtureBase.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/EdgeIterator.h>

#include <util/global.h>

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /*
   * Constructor.
   */
   Mixture::Mixture(MixtureBase const & mixture)
    : mixturePtr_(&mixture)
   {}

   /*
   * Destructor.
   */
   Mixture::~Mixture()
   {}

   /*
   * Allocate memory.
   */
   void Mixture::allocate()
   {
      const int nMonomer = mixture().nMonomer();
      const int nPolymer = mixture().nPolymer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nPolymer > 0);
      polymers_.allocate(nPolymer);
      for (int i = 0; i < nPolymer; ++i) {
         polymers_[i].associate(mixture().polymerSpecies(i));
         polymers_[i].allocate(nMonomer);
      }
   }

   /*
   * Compute mutable state variables.
   */
   void Mixture::setup()
   {
      const int nMonomer = mixture().nMonomer();
      const int nPolymer = mixture().nPolymer();

      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         kuhn[i] = mixture().monomer(i).kuhn();
      }

      for (int i = 0; i < nPolymer; ++i) {
         polymers_[i].setup(kuhn);
      }
   }

   /*
   * Compute k-space array for intramolecular correlation function.
   */
   void Mixture::computeOmega(int ma, int mb, 
                              Array<double> const & kSq, 
                              Array<double> & correlations) const
   {
      // Constants
      const double vMonomer = mixture().vMonomer();
      const int nMonomer = mixture().nMonomer();
      const int nPolymer = mixture().nPolymer();
      const int nSolvent = mixture().nSolvent();
      const int nk = kSq.capacity();

      // Preconditions
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nPolymer > 0);
      UTIL_CHECK(nSolvent >= 0);
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(correlations.capacity() == nk);

      // Initialize correlations to zero
      for (int i = 0; i < nk; ++i){
         correlations[i] = 0.0;
      }

      double phi, length, c;
      int nA, nB, ia, ib, ja, jb;

      // Loop over polymer species
      for (int ip = 0; ip < nPolymer; ip++){

         // Polymer properties
         Polymer const & polymer = polymers_[ip];
         phi = polymer.phi();
         length = polymer.totalLength();
         c = phi/(length*vMonomer);

         if (ma == mb) {
            // Monomers of the same type
            nA = polymer.blockIds(ma).size(); // # of blocks of type ma
            for (ia = 0; ia < nA; ++ia) {
                ja = polymer.blockIds(ma)[ia]; // block index
                polymer.computeOmega(ja, ja, c, kSq, correlations);
            }
            if (nA > 1) {
               // If the polymer contains multiple blocks of type ma,
               // compute interblock correlations
               for (ia = 0; ia < nA; ++ia) {
                  ja = polymer.blockIds(ma)[ia]; // block index
                  for (ib = 0; ib < ia; ++ib) {
                     jb = polymer.blockIds(mb)[ib]; // block index
                     polymer.computeOmega(ja, jb, 2.0*c, kSq, correlations);
                  }
               }
            }
         } else {
            // Monomers of the different types
            nA = polymer.blockIds(ma).size(); // # of blocks of type ma
            nB = polymer.blockIds(mb).size(); // # of blocks of type mb
            for (ia = 0; ia < nA; ++ia) {
               ja = polymer.blockIds(ma)[ia]; // block index, monomer type ma
               for (ib = 0; ib < nB; ++ib) {
                  jb = polymer.blockIds(mb)[ib]; // block index, monomer type mb
                  polymer.computeOmega(ja, jb, c, kSq, correlations);
               }
            }
         }

      } // loop over polymer species

      // Loop over solvent species (if any)
      if (ma == mb && nSolvent > 0) {
         double size;
         int monomerId;
         for (int i = 0; i < nSolvent; i++) {
            SolventSpecies const & solvent = mixture().solventSpecies(i);
            monomerId = solvent.monomerId();
            if (monomerId == ma) {
               phi = solvent.phi();
               size = solvent.size();
               c = phi*size/vMonomer;
               for (int i = 0; i < nk; ++i) {
                  correlations[i] += c;
               }
            }
         }
      }

   }

   /*
   * Compute k-space array of total intramolecular correlation function.
   */
   void Mixture::computeOmegaTotal(Array<double> const & kSq, 
                                   Array<double> & correlations) const
   {
      // Constants
      const double vMonomer = mixture().vMonomer();
      const int nPolymer = mixture().nPolymer();
      const int nSolvent = mixture().nSolvent();
      const int nk = kSq.capacity();

      // Preconditions
      UTIL_CHECK(nPolymer > 0);
      UTIL_CHECK(nSolvent >= 0);
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(correlations.capacity() == nk);

      // Initialize correlations to zero
      for (int i = 0; i < nk; ++i){
         correlations[i] = 0.0;
      }

      // Loop over polymer species
      double phi, length, c;
      for (int ip = 0; ip < nPolymer; ip++){
         Correlation::Polymer const & polymer = polymers_[ip];
         phi = polymer.phi();
         length = polymer.totalLength();
         c = phi/(length*vMonomer);
         polymer.computeOmegaTotal(c, kSq, correlations);
      } 

      // Loop over solvent species (if any)
      if (nSolvent > 0) {
         double size;
         for (int i = 0; i < nSolvent; i++){
            SolventSpecies const & solvent = mixture().solventSpecies(i);
            phi = solvent.phi();
            size = solvent.size();
            c = phi*size/vMonomer;
            for (int i = 0; i < nk; ++i) {
               correlations[i] += c;
            }
         }
      }

   }

}
}
#endif

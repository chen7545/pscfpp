#ifndef PSCF_OMEGA_MIXTURE_H
#define PSCF_OMEGA_MIXTURE_H

/*
* PSCF - Mixture Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/correlation/Polymer.h>
#include <util/containers/DArray.h>

// Forward declarations
namespace Pscf {
   class MixtureBase;
}

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /**
   * Intramolecular correlation analysis for one polymer Species
   *
   * \ingroup Pscf_Omega__Module
   */
   class Mixture 
   {

   public:

      /**
      * Constructor.
      *
      * \param mixture parent System object
      */
      Mixture(MixtureBase const& mixture);

      /**
      * Destructor.
      */
      ~Mixture();

      /*
      * Allocate private data structures, construct immutable data.
      *
      * The number of monomers and number of polymer species must be
      * known before this function is invoked.
      */
      void allocate();

      /*
      * Compute mutable private data.
      *
      * This function compute properties that may depend on block
      * lengths and statistical segment lengths.
      */
      void setup();

      /**
      * Compute ideal gas correlation function for a monomer type pair.
      */
      void computeOmega(int ma, int mb, 
                        Array<double> const & kSq,
                        Array<double> & correlation) const;

      /**
      * Compute total ideal gas density correlation function.
      */
      void computeOmegaTotal(Array<double> const & kSq,
                             Array<double> & correlation) const;

      /** 
      * Return reference to associated mixture descriptor.
      *
      * \param i  index for polymer species
      */      
      Correlation::Polymer const & polymer(int i) const;
      
   protected:

      /** 
      * Return reference to parent mixture.
      */      
      MixtureBase const & mixture() const;
      
   private:

      /// Array of Correlation::Polymer objects
      DArray<Correlation::Polymer> polymers_;
   
      /// Pointer to the associated MixtureBase object.
      MixtureBase const * mixturePtr_;
  
   };

   // Get a constituent Correlation::Polymer object by reference.
   inline Correlation::Polymer const & Mixture::polymer(int i) const
   {
      UTIL_CHECK(polymers_.isAllocated());
      return polymers_[i]; 
   }
   
   // Get the parent mixture by const reference.
   inline MixtureBase const & Mixture::mixture() const
   {  return *mixturePtr_; }
   
} // namespace Correlation
} // namespace Pscf
#endif

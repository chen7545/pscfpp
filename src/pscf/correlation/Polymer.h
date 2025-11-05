#ifndef PSCF_OMEGA_POLYMER_H
#define PSCF_OMEGA_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <util/containers/GArray.h>

// Forward declarations
namespace Pscf {
   class PolymerSpecies;
}

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /**
   * Intramolecular correlation analysis for one polymer Species
   */
   class Polymer 
   {

   public:

      /**
      * Default constructor.
      */
      Polymer();

      /**
      * Constructor (creates association).
      *
      * \param polymer  associated PolymerSpecies object
      */
      Polymer(PolymerSpecies const& polymer);

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Create association.
      *
      * \param polymer  associated PolymerSpecies object
      */
      void associate(PolymerSpecies const& polymer);

      /**
      * Allocate private memory.
      *
      * \pre nBlock must have been set for associated PolymerSpecies
      */
      void allocate(int nMonomer);

      /**
      * Setup private data.
      */
      void setup(Array<double> const & kuhn);

      /**
      * Compute intramolecular correlation function for a pair of blocks.
      *
      * This function evaluates the intramolecular correlation function
      * Omega for a pair of blocks with block indices ia and ib at a 
      * value of the squared wavenumber given by input parameter kSq, 
      * and returns the result. The return value is given by the product 
      * of a dimensionless single-chain correlation function omega and the 
      * input parameter prefactor. To obtain the contribution of this 
      * species to the correlation function for a homogeneous ideal gas
      * mixture, set prefactor equal to the species chain concentration 
      * prefactor = phi/(v*totalLength), where v is the monomer reference 
      * volume. 
      *
      * \param ia  block index of first block
      * \param ib  block index of second block
      * \param prefactor  prefactor multiplying omega(k)
      * \param kSq  squared wavenumber 
      * \return calculated value of Omega(k)
      */
      double computeOmega(int ia, int ib, double prefactor, double kSq) 
      const;

      /**
      * Compute intramolecular correlation function for a pair of blocks.
      * 
      * This function evaluates the intramolecular correlation function
      * Omega for a pair of blocks with block indices ia and ib at a list
      * of values for the squared wavenumber give by elements of input
      * array parameter kSq and adds each result to the corresponding 
      * element of array correlation. The value for each wavenumber is 
      * given by the product of a dimensionless single-chain correlation 
      * function omega(k) and the input parameter "prefactor". To obtain 
      * the contribution of this species to correlation function for an
      * ideal gas mixture, set prefactor = phi/(v*totalLength), where v 
      * is the monomer reference volume. 
      *
      * Resulting values of Omega are added to elements of the array 
      * correlation. This is to calculations in which results of multiple 
      * polymer species are added. All elements of array correlation 
      * should thus be set to zero before beginning such a calculation.
      *
      * \param ia  block index of first block
      * \param ib  block index of second block
      * \param prefactor  prefactor multiplying omega(q)
      * \param kSq  array of squared wavenumbers (in)
      * \param correlation  array of correlation functions (out)
      */
      void computeOmega(int ia, int ib, double prefactor,
                        Array<double> const & kSq, 
                        Array<double> & correlation) const;

      /**
      * Compute total intramolecular correlation function.
      * 
      * This function computes the total density-density correlation 
      * function for a list of values for the squared wavenumber give
      * by elements of input array parameter kSq and adds the result 
      * for each to the corresponding element of array "correlation". 
      * The value for each wavenumber is given by the product of a 
      * dimensionless single-chain correlation function omega(k) and
      * the "prefactor" input parameter. To obtain the contribution 
      * of this species to the correlation function for an ideal gas 
      * mixture, set the prefactor = phi/(v*totalLength), where v is 
      * the monomer reference volume. 
      *
      * Resulting values of Omega are added to elements of output array 
      * correlation.
      *
      * \param prefactor  prefactor multiplying omega(q)
      * \param kSq  array of squared wavenumbers (in)
      * \param correlation  array of correlation functions (out)
      */
      void computeOmegaTotal(double prefactor,
                             Array<double> const & kSq, 
                             Array<double> & correlation) const;

      /**
      * Return the volume fraction of this polymer species.
      */
      double phi() const;

      /**
      * Return the sum of lengths of all blocks in this polymer species.
      */
      double totalLength() const;

      /**
      * Get the number of blocks in the associated polymer species.
      */
      int nBlock() const;

      /**
      * Get a list of block ids for blocks of specified monomer type.
      *
      * \param i  monomer type index
      */
      GArray<int> const & blockIds(int i) const;

      /**
      * Get the mean-squared length of path between blocks i and j.
      *
      * \param i  block index of 1st block.
      * \param j  block index of 2nd block.
      */
      double rSq(int i, int j) const;

   private:
     
      /*
      * Block lengths, indexed by block index.
      */
      DArray<double> length_;

      /*
      * Monomer statistical segment lengths, indexed by block index.
      */
      DArray<double> kuhn_;
    
      /* 
      * Symmetric matrix of values of rSq for paths between blocks.
      */
      DMatrix<double> rSq_;

      /*
      * Identifiers for blocks of specified monomer type.
      *
      * Elements of the outer DArray are indexed by monomer type id. 
      * Each element is a GArray of blocks ids for all blocks of one 
      * monomer type.
      */
      DArray< GArray<int> > blockIds_;

      /* 
      * Pointer to the associated PolymerSpecies object.
      */
      PolymerSpecies const * speciesPtr_;

      /// Volume fraction of this polymer
      double phi_;

      /// Sum of lengths of all blocks in the polymer
      double totalLength_;

      /// Number of blocks in this polymer.
      int nBlock_;

      /// Number of monomer types in the mixture.
      int nMonomer_;

      /** 
      * Return reference to associated polymer species.
      */      
      PolymerSpecies const & species() const;
 
   };
  
   // Public inline member functions
 
   // Get the volume fraction of this polymer species.
   inline double Polymer::phi() const
   {  return phi_; }
   
   // Get the sum of the lengths of all blocks in this polymer species.
   inline double Polymer::totalLength() const
   {  return totalLength_; }
   
   // Get the number of blocks in this polymer species.
   inline int Polymer::nBlock() const
   {  return nBlock_; }
   
   // Get the list of blocks ids for a specified monomer type.
   inline GArray<int> const & Polymer::blockIds(int i) const
   {  return blockIds_[i]; }
   
   // Get the mean-squared length of the path connecting two blocks.
   inline double Polymer::rSq(int i, int j) const
   {  return rSq_(i, j); }
   
   // Private inline member function
 
   // Get the associated PolymerSpecies molecule descriptor object.
   inline PolymerSpecies const & Polymer::species() const
   {  return *speciesPtr_; }
   
} // namespace Correlation
} // namespace Pscf
#endif

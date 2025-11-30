#ifndef CPC_INTERACTION_H
#define CPC_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/Matrix.h>       // argument (template)

namespace Pscf {
namespace Cpc {

   using namespace Util;

   /**
   * Interaction model for complex Langevin FTS.
   *
   * \ingroup Cpc_System_Module
   */
   class Interaction : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Interaction();

      /**
      * Destructor.
      */
      virtual ~Interaction();

      // Modifiers

      /**
      * Set the number of monomer types.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Read model parameters.
      *
      * \pre Must be called after setNMonomer.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Change one element of the chi matrix.
      *
      * \param i row index
      * \param j column index
      * \param chi  input value of chi
      */
      void setChi(int i, int j, double chi);

      // Accessors

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      /**
      * Return the chi matrix by const reference.
      */
      DMatrix<double> const & chi() const;

      /**
      * Return one element of the chi matrix.
      *
      * \param i row index
      * \param j column index
      */
      double chi(int i, int j) const;

      /** 
      * Return the dimensionless compression modulus.
      */
      double zeta() const; 

      /** 
      * Return the range of binary interactions.
      */
      double range() const; 

      /**
      * Compute the Fourier-space damping factor for pair interactions.
      *
      * \param kSq  squared wavenumber (input)
      * \return damping Fourier-space factor for pair interactions
      */
      double g(double kSq) const;

   private:

      // Symmetric matrix of interaction parameters.
      DMatrix<double> chi_;

      // Dimensionless compression modulus (a la Helfand)
      double zeta_;

      // Range of interaction
      double range_;

      // Constant used in computation of g(k)
      double alpha_;

      // Number of monomers
      int nMonomer_;

   };

   // Inline function

   inline int Interaction::nMonomer() const
   {  return nMonomer_; }

   inline DMatrix<double> const &  Interaction::chi() const
   {  return chi_; }

   inline double Interaction::chi(int i, int j) const
   {  return chi_(i, j); }

   inline double Interaction::zeta() const
   {  return zeta_; }

   inline double Interaction::range() const
   {  return range_; }

} // namespace Cpc
} // namespace Pscf
#endif

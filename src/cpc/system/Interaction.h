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
      * 
      * If the isCompressible parameter of this contructor is false, the 
      * physical model will be assumed to be incompressible, in which case 
      * a zeta parameter may not appear in the readParameters function. 
      * If this constructor parameter is true, then zeta is treated as 
      * an optional parameter in the parameter file, and the model is 
      * taken to be compressible if a value for zeta is given in the 
      * parameter file, and incompressible otherwise.
      *
      * \param isCompressible  allow use of a finite zeta parameter
      */
      Interaction(bool isCompressible = false);

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
      * Is the system compressible?
      *
      * The boolean flag that is returned by this function is initially
      * set equal to the value of the parameter passed to the constructor.
      * If this initial value is false, it may not be changed later. If
      * this initial value is true, then the value it may be reset in the
      * readParameters method, which sets isCompressible false if the
      * parameter file does not contain a zeta parameter.
      */
      bool isCompressible() const; 

   private:

      // Symmetric matrix of interaction parameters.
      DMatrix<double> chi_;

      // Dimensionless compression modulus (a la Helfand)
      double zeta_;

      // Number of monomers.
      int nMonomer_;

      // True if a finite value of zeta is defined.
      bool isCompressible_;

      // Have the object been initialized via readParameters.
      bool isInitialized_;

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

   inline bool Interaction::isCompressible() const
   {  return isCompressible_; }

} // namespace Cpc
} // namespace Pscf
#endif

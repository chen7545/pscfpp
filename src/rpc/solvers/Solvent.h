#ifndef RPC_SOLVENT_H
#define RPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>   // base class
#include <prdc/cpu/RField.h>            // member

// Forward declaration
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref user_param_solvent_sec "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Solvent : public SolventSpecies<double>
   {
   public:

      /**
      * Constructor.
      */
      Solvent();

      /**
      * Destructor.
      */
      ~Solvent();

      /**
      * Create an association with a mesh.
      *
      * Must be called before allocate().
      *
      * \param mesh  associated Mesh<D> object - spatial discretization
      */
      void associate(Mesh<D> const & mesh);

      /**
      * Allocate memory for a concentration field.
      *
      * This function may only be called once during initialization,
      * after the associate() function and before the first call to
      * the compute() function.
      */
      void allocate();

      /**
      * Compute concentration field, q, and phi or mu.
      *
      * Computes monomer concentration field cField, partition function
      * q, and either the solvent volume fraction phi or solvent chemical
      * potential mu, depending on ensemble. The function takes the
      * chemical potential field wField for the relevant monomer type as
      * its only input argument.
      *
      * The optional parameter phiTot is only relevant to problems such
      * as thin films in which the material is excluded from part of the
      * unit cell by imposing an inhomogeneous constraint on the sum of
      * monomer concentrations (i.e., a "mask").
      *
      * \param wField  monomer chemical potential field of relevant type.
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(RField<D> const & wField, double phiTot = 1.0);

      /**
      * Get the monomer concentration field for this solvent.
      */
      RField<D> const & cField() const;

   private:

      /// Concentration field for this solvent.
      RField<D> cField_;

      /// Pointer to associated mesh.
      Mesh<D> const *  meshPtr_;

   };

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D> inline
   RField<D> const & Solvent<D>::cField() const
   {  return cField_;  }

   // Explicit instantiation declarations
   extern template class Solvent<1>;
   extern template class Solvent<2>;
   extern template class Solvent<3>;

}
}
#endif

#ifndef CPC_MIXTURE_H
#define CPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cl/Mixture.h>     // base class template
#include <cpc/system/Types.h>    // base class argument

namespace Pscf {
namespace Cpc {

   // Forward declarations
   template <int D> class Polymer;
   template <int D> class Solvent;

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from a partial specialization of the template
   * Prdc::Cl::Mixture, and has the same public interface as this base
   * class template.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Cpc_Solver_Module
   */
   template <int D>
   class Mixture 
    : public Cl::Mixture<D, Polymer<D>, Solvent<D>, Types<D> >
   {

   public:

      /// Direct base class, instantiation of Cl::Mixture template.
      using ClMixtureT
         = typename Prdc::Cl::Mixture<D, Polymer<D>, Solvent<D>, Types<D> >;

      // Inherited public type name aliases

      using typename ClMixtureT::MixtureTmplT;
      using typename ClMixtureT::MixtureBaseT;
      using typename ClMixtureT::PolymerT;
      using typename ClMixtureT::SolventT;
      using typename ClMixtureT::BlockT;
      using typename ClMixtureT::PropagatorT;
      using typename ClMixtureT::FieldT;
      using typename ClMixtureT::FFTT;
      using typename ClMixtureT::WaveListT;

      // Inherited public member functions

      using ClMixtureT::readParameters;
      using ClMixtureT::associate;
      using ClMixtureT::allocate;
      using ClMixtureT::clearUnitCellData;
      using ClMixtureT::setKuhn;
      using ClMixtureT::compute;

      using MixtureTmplT::polymer;
      using MixtureTmplT::polymerSpecies;
      using MixtureTmplT::solvent;
      using MixtureTmplT::solventSpecies;

      using MixtureBaseT::nMonomer;
      using MixtureBaseT::monomer;
      using MixtureBaseT::nPolymer;
      using MixtureBaseT::nSolvent;
      using MixtureBaseT::nBlock;
      using MixtureBaseT::vMonomer;
      using MixtureBaseT::isCanonical;

   protected:

      using ClMixtureT::mesh;
      using ClMixtureT::ds;

   private:

      /**
      * Set all elements of a field to a common scalar: A[i] = s.
      *
      * \param A  field (LHS)
      * \param s  scalar (RHS)
      */
      void eqS(FieldT& A, double s) const override;

      /**
      * Compound addition assignment for fields : A[i] += B[i].
      *
      * \param A  field (LHS)
      * \param B  field (RHS)
      */
      void addEqV(FieldT& A, FieldT const & B) const override;

      /**
      * Allocate memory for all blocks
      */
      void allocateBlocks() override;

   };

   // Explicit instantiation declarations for derived class
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;

} // namespace Cpc
namespace Prdc {

   // Explicit instantiation declarations for base class
   extern template 
   class Cl::Mixture<1, Cpc::Polymer<1>, Cpc::Solvent<1>, Cpc::Types<1> >;
   extern template 
   class Cl::Mixture<2, Cpc::Polymer<2>, Cpc::Solvent<2>, Cpc::Types<2> >;
   extern template 
   class Cl::Mixture<3, Cpc::Polymer<3>, Cpc::Solvent<3>, Cpc::Types<3> >;

} // namespace Prdc
} // namespace Pscf
#endif

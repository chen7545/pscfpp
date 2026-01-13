#ifndef RPC_MIXTURE_H
#define RPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Mixture.h>   // base class template
#include <rpc/system/Types.h>     // base class template argument

namespace Pscf {
namespace Rpc {

   // Forward declarations
   template <int D> class Polymer;
   template <int D> class Solvent;

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from an instantiation of the class template
   * Rp::Mixture, and has the same public interface as this base class
   * class template.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Mixture : public Rp::Mixture<D, Polymer<D>, Solvent<D>, Types<D> >
   {

   public:

      /// Direct (parent) base class.
      using RpMixtureT
         = typename Rp::Mixture<D, Polymer<D>, Solvent<D>, Types<D> >;

      // Inherited public type name aliases

      using typename RpMixtureT::MixtureTmplT;
      using typename RpMixtureT::PolymerT;
      using typename RpMixtureT::SolventT;
      using typename RpMixtureT::BlockT;
      using typename RpMixtureT::PropagatorT;
      using typename RpMixtureT::FieldT;
      using typename RpMixtureT::FFTT;
      using typename RpMixtureT::WaveListT;

      // Inherited public member functions

      using RpMixtureT::readParameters;
      using RpMixtureT::associate;
      using RpMixtureT::allocate;
      using RpMixtureT::clearUnitCellData;
      using RpMixtureT::setKuhn;
      using RpMixtureT::compute;
      using RpMixtureT::computeStress;
      using RpMixtureT::hasStress;
      using RpMixtureT::createBlockCRGrid;

      using MixtureTmplT::polymer;
      using MixtureTmplT::polymerSpecies;
      using MixtureTmplT::solvent;
      using MixtureTmplT::solventSpecies;

      using MixtureBase<double>::nMonomer;
      using MixtureBase<double>::monomer;
      using MixtureBase<double>::nPolymer;
      using MixtureBase<double>::nSolvent;
      using MixtureBase<double>::nBlock;
      using MixtureBase<double>::vMonomer;
      using MixtureBase<double>::isCanonical;

   protected:

      using RpMixtureT::mesh;
      using RpMixtureT::ds;

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

} // namespace Rpc

namespace Rp {
   // Explicit instantiation declarations for base class
   extern template 
   class Mixture<1, Rpc::Polymer<1>, Rpc::Solvent<1>, Rpc::Types<1> >;
   extern template 
   class Mixture<2, Rpc::Polymer<2>, Rpc::Solvent<2>, Rpc::Types<2> >;
   extern template 
   class Mixture<3, Rpc::Polymer<3>, Rpc::Solvent<3>, Rpc::Types<3> >;
} 

} // namespace Pscf
#endif

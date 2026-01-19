#ifndef RPG_MIXTURE_H
#define RPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Mixture.h>    // base class template
#include <rpg/system/Types.h>      // base class template argument
#include <pscf/cuda/cudaTypes.h>   // real & complex cuda data types

namespace Pscf {
namespace Rpg {

   // Forward declarations
   template <int D> class Polymer;
   template <int D> class Solvent;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from an instantiation of the class template
   * Rp::Mixture, and has the same public interface as this base class
   * template. 
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Mixture : public Rp::Mixture<D, Types<D> >
   {

   public:

      /// Direct (parent) base class.
      using RpMixtureT = typename Rp::Mixture<D, Types<D> >;

      // Inherited public type name aliases

      using typename RpMixtureT::MixtureTmplT;
      using typename RpMixtureT::PolymerT;
      using typename RpMixtureT::SolventT;
      using typename RpMixtureT::BlockT;
      using typename RpMixtureT::PropagatorT;
      using typename RpMixtureT::FieldT;
      using typename RpMixtureT::FFTT;
      using typename RpMixtureT::WaveListT;

      // Public member functions

      /**
      * Constructor.
      */
      Mixture();

      /**
      * Read all parameters and initialize.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in) override;

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

      using MixtureBase<cudaReal>::nMonomer;
      using MixtureBase<cudaReal>::monomer;
      using MixtureBase<cudaReal>::nPolymer;
      using MixtureBase<cudaReal>::nSolvent;
      using MixtureBase<cudaReal>::nBlock;
      using MixtureBase<cudaReal>::vMonomer;
      using MixtureBase<cudaReal>::isCanonical;

   protected:

      using RpMixtureT::mesh;
      using RpMixtureT::ds;

   private:

      // Private member data

      /// Use batched FFTs to compute stress? (faster, but doubles memory)
      bool useBatchedFFT_;

      // Private member functions

      /**
      * Set all elements of a field to a common scalar: A[i] = s.
      *
      * \param A  field (LHS)
      * \param s  scalar (RHS)
      */
      virtual void eqS(FieldT& A, double s) const override;

      /**
      * Compound addition assignment for fields : A[i] += B[i].
      *
      * \param A  field (LHS)
      * \param B  field (RHS)
      */
      virtual void addEqV(FieldT& A, FieldT const & B) const override;

      /**
      * Allocate memory for all blocks
      */
      virtual void allocateBlocks() override;

   };

   // Explicit instantiation declarations
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;

} // namespace Rpg

namespace Rp {
   // Explicit instantiation declarations for base class
   extern template class Mixture<1, Rpg::Types<1> >;
   extern template class Mixture<2, Rpg::Types<2> >;
   extern template class Mixture<3, Rpg::Types<3> >;

} // namespace Prdc

} // namespace Pscf
#endif

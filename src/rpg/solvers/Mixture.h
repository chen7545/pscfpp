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

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from an instantiation of the class template
   * Rp::Mixture, and has the same public interface as this base class
   * template. See documentation of the base class for details.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Mixture : public Rp::Mixture< D, Types<D> >
   {

   public:

      /// Direct (parent) base class.
      using RpMixtureT = typename Rp::Mixture<D, Types<D> >;

      // Inherited public names (for use in member functions)
      using typename RpMixtureT::MixtureTmplT;
      using typename RpMixtureT::MixtureBaseT;
      using typename RpMixtureT::FieldT;
      using MixtureTmplT::polymer;

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

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations for derived and base classes
namespace Pscf {
   extern template class MixtureTmpl< Rpg::Polymer<1>, Rpg::Solvent<1> >;
   extern template class MixtureTmpl< Rpg::Polymer<2>, Rpg::Solvent<2> >;
   extern template class MixtureTmpl< Rpg::Polymer<3>, Rpg::Solvent<3> >;
   namespace Rp {
      extern template class Mixture<1, Rpg::Types<1> >;
      extern template class Mixture<2, Rpg::Types<2> >;
      extern template class Mixture<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Mixture<1>;
      extern template class Mixture<2>;
      extern template class Mixture<3>;
   }
}
#endif

#ifndef RPG_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RPG_EINSTEIN_CRYSTAL_PERTURBATION_H

#include "Perturbation.h"            // base class
#include <prdc/cuda/RField.h>        // member
#include <util/containers/DArray.h>  // member

namespace Pscf {
namespace Rpg {

   template <int D> class Simulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Perturbation for Einstein crystal thermodynamic integration.
   *
   * \see \ref rp_EinsteinCrystalPerturbation_page "Einstein Crystal"
   * \see \ref psfts_perturb_page "Perturbations"
   * \ingroup Rpg_Fts_Perturbation_Module
   */
   template <int D>
   class EinsteinCrystalPerturbation : public Perturbation<D>
   {

   public:

      /**
      * Constructor.
      *  
      * \param Simulator  parent Simulator object
      */
      EinsteinCrystalPerturbation(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~EinsteinCrystalPerturbation();

      /**
      * Read body of parameter file block and initialize.
      *
      * \param in input parameter file stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      */
      virtual void setup();

      /**
      * Compute and return the perturbation to the Hamiltonian.
      *
      * \param unperturbedHamiltonian Hamiltonian without perturbation
      */
      virtual double hamiltonian(double unperturbedHamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      *
      * \param dc  functional derivatives of Hamiltonian (in/out)
      */
      virtual void incrementDc(DArray< RField<D> >& dc);

      /**
      * Compute and return derivative of free energy w/respect to lambda.
      */
      virtual double df();

      /**
      * Save any required internal state variables.
      */
      virtual void saveState();

      /**
      * Restore any required internal state variables.
      */
      virtual void restoreState();

   protected:

      // Inherited protected members
      using Perturbation<D>::lambda_;
      using Perturbation<D>::simulator;
      using Perturbation<D>::system;
      using ParamComposite::readOptionalDArray;

   private:

      // Parameters used in Einstein crystal integration
      DArray<double> epsilon_;

      // Reference w field
      DArray< RField<D> > w0_;

      // Eigenvector components of the reference w fields
      DArray< RField<D> > wc0_;

      // Work space
      RField<D> dw_;

      // Current Einstein crystal Hamiltonian
      double ecHamiltonian_;

      // Current unperturbed Hamiltonian
      double unperturbedHamiltonian_;

      // Saved Einstein crystal Hamiltonian
      double stateEcHamiltonian_;

      // Saved unperturbed Hamiltonian
      double stateUnperturbedHamiltonian_;

      // Have epsilon values been set?
      bool hasEpsilon_;

      // Reference field file name
      std::string fieldFileName_;

      // Compute eigenvector components of the reference field
      void computeWcReference();

   };

   // Explicit instantiation declarations
   extern template class EinsteinCrystalPerturbation<1>;
   extern template class EinsteinCrystalPerturbation<2>;
   extern template class EinsteinCrystalPerturbation<3>;

}
}
#endif

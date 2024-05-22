#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_H

#include "Perturbation.h"      // base class

namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Prdc{
namespace Cpu{
   template <int D> class RField;
}
}
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   template <int D> class Simulator;

   /**
   * Perturbation for Einstein crystal thermodynamic integration method.
   *
   * \ingroup Rpc_Simulate_Perturbation_Module
   */
   template <int D>
   class EinsteinCrystalPerturbation : public Perturbation<D>
   {

   public:

      /**
      * Constructor.
      */
      EinsteinCrystalPerturbation(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~EinsteinCrystalPerturbation();

      /**
      * Read parameters from archive.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      */
      virtual void setup();

      /**
      * Compute and return the perturbation to the Hamiltonian.
      */
      virtual double modifyHamiltonian(double hamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      */
      virtual void modifyDc(DArray< RField<D> >& dc);
      
      /**
      * Compute and return derivative of free energy.
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
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      
   protected:
      
      // Inherited protected functions
      using Perturbation<D>::simulator;
      using Perturbation<D>::system;
      using Perturbation<D>::lambda;
      using Perturbation<D>::mode;
      
      // Inherited protected data members
      using Perturbation<D>::lambda_;
      using Perturbation<D>::mode_;
      
    private:

      // Initial coupling parameter
      double lambda0_;
      
      // Increment rate
      double dLambda_;
      
      // Spring constant for the einstein crystal.
      double alpha_;
      
      // Reference w field
      DArray< RField<D> > w0_;
      
      // Eigenvector components of the reference w fields
      DArray< RField<D> > wc0_;
      
      // Current Einstein Crystal hamiltonian 
      double hamiltonianEC_;
      
      // Current Block copolymer hamiltonian 
      double hamiltonianBCP_;
      
      // Saved Einstein Crystal hamiltonian  
      double stateHamiltonianEC_;
      
      // Saved Block copolymer hamiltonian 
      double stateHamiltonianBCP_;
      
      // Reference FieldFileName
      std::string referenceFieldFileName_;
      
      // Compute eigenvector components of the reference field
      void computeWcReference();
   
   };

   #ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP
   // Suppress implicit instantiation
   extern template class EinsteinCrystalPerturbation<1>;
   extern template class EinsteinCrystalPerturbation<2>;
   extern template class EinsteinCrystalPerturbation<3>;
   #endif

}
}
#endif

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

   template <int D> class System;
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
      * Compute and return the perturbated hamiltonian.
      */
      virtual double modifyHamiltonian(double hamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      */
      virtual void modifyDc(DArray< RField<D> >& dc);
      
      void updateDf();
      
      /**
      * Compute Integral pver lambda
      */
      void computeIntegral();
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using Perturbation<D>::simulator;
      using Perturbation<D>::system;

   private:

      // Coupling parameter
      double lambda_;
      
      // Ramping rate
      double dLambda_;
      
      // Free energy 
      double f_;
      
      // derivative of free energy 
      double df_;
      
      // output interval 
      long interval_;
      
      // Reference w field
      DArray< RField<D> > w0_;
      
      // Eigenvector components of the reference w fields
      DArray< RField<D> > wc0_;
      
      // Einstein Crystal hamiltonian 
      double hamiltonianEC_;
      
      // Block copolymer hamiltonian 
      double hamiltonianBCP_;
      
      // Reference FieldFileName
      std::string referenceFieldFileName_;
      
      // Output file Name
      std::string outputFileName_;
      
      // Output file stream
      std::ofstream outputFile_;
      
      // 
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

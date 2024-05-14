#ifndef RPC_PERTURBATION_TPP
#define RPC_PERTURBATION_TPP

#include "Perturbation.h"
#include <rpc/simulate/Simulator.h>
#include <prdc/cpu/RField.h>

#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   Perturbation<D>::Perturbation(Simulator<D>& simulator)
    : ParamComposite(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {}
   
   /* 
   * Destructor.
   */
   template <int D>
   Perturbation<D>::~Perturbation()
   {}
   
   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::readParameters(std::istream& in)
   {}

   /*
   * Setup before simulation, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::setup()
   {}

   /*
   * Compute and return perturbation, default implementation.
   */
   template <int D>
   double Perturbation<D>::modifyHamiltonian(double hamiltonian)
   { return hamiltonian; }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::modifyDc(DArray< RField<D> > & dc)
   {}
   
   /*
   * Modify free energy derivatives, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::updateDf()
   {}
   
   /*
   * Return BCP Hamiltonian
   */
   template <int D>
   double Perturbation<D>::bcpHamiltonian()
   { return 0;}
   
   /*
   * Return BCP Hamiltonian
   */
   template <int D>
   double Perturbation<D>::ecHamiltonian()
   { return 0;}
   
   /**
   * Modify BCP Hamiltonian
   */
   template <int D>
   void Perturbation<D>::restoreBCPHamiltonian(double hamiltonian)
   {}
   
   /**
   * Modify EC Hamiltonian
   */
   template <int D>
   void Perturbation<D>::restoreECHamiltonian(double hamiltonian)
   {}
      
}
}
#endif 

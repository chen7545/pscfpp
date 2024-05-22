#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <util/param/ParamComposite.h>      // base class

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
   * Abstract base for perturbations of standard Hamiltonian.
   *
   * \ingroup Rpc_Simulate_Perturbation_Module
   */
   template <int D>
   class Perturbation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Perturbation(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~Perturbation();

      /**
      * Read parameters from archive.
      *
      * Empty default implementation.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      *
      * This method must be called just before the beginning of
      * the main simulation loop, after an initial configuration 
      * is known. It may be used to complete any initialization
      * that cannot be completed in the readParam method, because
      * knowledge of the configuration is needed. 
      *
      * The default implementation is an empty function.
      */
      virtual void setup();

      /**
      * Modify and return Hamiltonian to include perturbation.
      *
      * Default implementation returns unperturbated hamiltonian. 
      */
      virtual double modifyHamiltonian(double hamiltonian);
      
      /**
      * Modify the generalized forces to include perturbation.
      *
      * Empty default implementation.
      */
      virtual void modifyDc(DArray< RField<D> >& dc);
      
      /**
      * Compute and return derivative of free energy 
      * with respective to specific variable per monomer.
      * 
      * Default implementation returns 0. 
      */ 
      virtual double df();

      /**
      * Save any required internal state variables.
      *
      * This function should save any state variables that would need to 
      * be restored after a rejected Monte Carlo move or failure of the
      * compressor to converge after an attempted Brownian dynamics move.
      *
      * Empty default implementation 
      */
      virtual void saveState()
      {};

      /**
      * Restore any required internal state variables.
      *
      * This function is called after rejection of an MC move or failure
      * of an attempted BD step, and should restore the variables saved 
      * by the saveState function.
      *
      * Empty default implementation 
      */
      virtual void restoreState()
      {};
      
      /**
      * Get the strength of the perturbation
      */
      double lambda() const;
      
      /**
      * Get the mode of the thermodynamic integration
      */
      int mode() const;

      /**
      * Get parent Simulator<D> by const reference.
      */
      Simulator<D> const & simulator() const;
      
      /** 
      * Get parent System<D> by non-const reference.
      */      
      System<D>& system();

   protected:
      
      /**
      * Strength of the perturbation
      */
      double lambda_;
   
      /**
      * The mode of thermodynamic integration
      * mode = 0 correspond to a static parameter
      * mode = 1 be continuous.
      * 
      */
      int mode_;
      
      /**
      * Get parent Simulator<D> by non-const reference.
      */
      Simulator<D>& simulator();

   private:

      /// Pointer to parent Simulator.
      Simulator<D>* simulatorPtr_;
      
      /// Pointer to parent System.
      System<D>* systemPtr_;

   };

   // Inline methods
   
   // Get the strength of the perturbation
   template <int D>
   inline double Perturbation<D>::lambda() const
   { return lambda_; }
   
   // Get the mode of the thermodynamic integration
   template <int D>
   inline int Perturbation<D>::mode() const
   { return mode_; }

   // Return parent simulator by const reference.
   template <int D>
   inline Simulator<D> const & Perturbation<D>::simulator() const
   {
      assert(simulatorPtr_);  
      return *simulatorPtr_; 
   }

   // Return parent simulator by non-const reference.
   template <int D>
   inline Simulator<D> & Perturbation<D>::simulator() 
   {  
      assert(simulatorPtr_);  
      return *simulatorPtr_; 
   }
   
   // Return parent simulator by non-const reference.
   template <int D>
   inline System<D> & Perturbation<D>::system() 
   {  
      assert(systemPtr_);  
      return *systemPtr_; 
   }

   // Method template

   #ifndef RPC_PERTURBATION_TPP
   // Suppress implicit instantiation
   extern template class Perturbation<1>;
   extern template class Perturbation<2>;
   extern template class Perturbation<3>;
   #endif

}
}
#endif

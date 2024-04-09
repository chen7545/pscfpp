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
      * Compute and return the perturbated to the Hamiltonian.
      *
      * Default implementation returns . 
      */
      virtual double modifyHamiltonian(double hamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      *
      * Empty default implementation.
      */
      virtual void modifyDc(DArray< RField<D> >& dc);
      
      virtual void updateDf();

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

   #ifndef PC_PERTURBATION_TPP
   // Suppress implicit instantiation
   extern template class Perturbation<1>;
   extern template class Perturbation<2>;
   extern template class Perturbation<3>;
   #endif

}
}
#endif

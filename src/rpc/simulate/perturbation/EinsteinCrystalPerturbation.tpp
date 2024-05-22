#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP

#include "EinsteinCrystalPerturbation.h"
#include <rpc/simulate/Simulator.h>
#include <rpc/System.h>
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
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(Simulator<D>& simulator)
    : Perturbation<D>(simulator),
      lambda0_(0.0),
      dLambda_(0.0),
      alpha_(0.0),
      hamiltonianEC_(0.0),
      hamiltonianBCP_(0.0),
      stateHamiltonianEC_(0.0),
      stateHamiltonianBCP_(0.0)
   { setClassName("EinsteinCrystal"); }
   
   /* 
   * Destructor.
   */
   template <int D>
   EinsteinCrystalPerturbation<D>::~EinsteinCrystalPerturbation()
   {}
   
   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::readParameters(std::istream& in)
   {
      // Readin
      read(in, "referenceFieldFileName", referenceFieldFileName_);
      read(in, "lambda", lambda0_);
      read(in, "alpha", alpha_);
      readOptional(in, "incrementRate", dLambda_);
      readOptional(in, "mode", mode_);
   }

   /*
   * Setup before simulation.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::setup()
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const IntVec<D> dimensions = system().domain().mesh().dimensions();
      
      UTIL_CHECK(nMonomer == 2);
      
      // Check mode and increment rate
      if (mode_ == 0){
         UTIL_CHECK(dLambda_ == 0.0);
      } else if (mode_ == 1){
         UTIL_CHECK(dLambda_ != 0.0);
      } else {
         UTIL_THROW("Invalid thermodynamic integration mode type.");
      }
      
      // Allocate memory for reference field
      w0_.allocate(nMonomer);
      wc0_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0_[i].allocate(dimensions);
         wc0_[i].allocate(meshSize);
      }
      
      // Read in reference field 
      system().fieldIo().readFieldsRGrid(referenceFieldFileName_, w0_, 
                                         system().domain().unitCell());
      
      // Compute eigenvector components of the reference field
      computeWcReference();
      
   }

   /*
   * Compute and return perturbation to Hamiltonian.
   */
   template <int D>
   double EinsteinCrystalPerturbation<D>::modifyHamiltonian(double hamiltonian)
   {
      lambda_ = lambda0_ + simulator().iStep()* dLambda_;
      
      // Compute EC Hamiltonian
      double prefactor, w, s;
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      hamiltonianEC_ = 0.0;
      
      for (int j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & Wc = simulator().wc(j);
         prefactor = alpha_;
         s = simulator().sc(j);
         for (int i = 0; i < meshSize; ++i) {
            w = Wc[i] - s - wc0_[j][i];
            hamiltonianEC_ += prefactor*w*w;
         }
      }
      
      // Normalize hamiltonianEC to equal a value per monomer
      hamiltonianEC_/= double(meshSize);
      hamiltonianEC_ *= nMonomerSystem;
      
      // Obtain hamiltonianBCP
      hamiltonianBCP_ = hamiltonian;  
  
      return lambda_* hamiltonianBCP_ + (1.0-lambda_) * hamiltonianEC_;
   }
   
   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void 
   EinsteinCrystalPerturbation<D>::modifyDc(DArray< RField<D> > & dc)
   {
      double DcBCP, DcEC;
      double prefactor, s;
      const int meshSize = system().domain().mesh().size();
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      
      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         RField<D>& Dc = dc[i];
         RField<D> const & Wc = simulator().wc(i);
         
         // Loop over grid points
         for (int k = 0; k < meshSize; ++k) {
            // Copy block copolymer derivative
            DcBCP = Dc[k];
            
            // Compute EC derivative
            prefactor = -1.0*double(nMonomer)/simulator().chiEval(i)/vMonomer;
            s = simulator().sc(i);
            DcEC = prefactor * (Wc[k] - s - wc0_[i][k]);
            
            // Compute composite derivative
            Dc[k] = lambda_* DcBCP + (1.0 - lambda_) * DcEC;
         }
      }
      
   }
   
   /**
   * Compute and return derivative of free energy.
   */ 
   template <int D>
   double EinsteinCrystalPerturbation<D>::df(){
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      
      if (mode_ == 1){
         return (hamiltonianBCP_ - hamiltonianEC_) * dLambda_/nMonomerSystem;
      }
      
      return (hamiltonianBCP_ - hamiltonianEC_)/nMonomerSystem;
   }
      
   /**
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::saveState(){
      stateHamiltonianEC_ = hamiltonianEC_;
      stateHamiltonianBCP_ = hamiltonianBCP_;
   }
   
   /**
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::restoreState(){
      hamiltonianEC_ = stateHamiltonianEC_;
      hamiltonianBCP_ = stateHamiltonianBCP_;
   }
   
   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs of the projected chi matrix as a basis.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::computeWcReference()
   {

      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Loop over eigenvectors (j is an eigenvector index)
      for (j = 0; j < nMonomer; ++j) {

         // Loop over grid points to zero out field wc_[j]
         RField<D>& Wc = wc0_[j];
         for (i = 0; i < meshSize; ++i) {
            Wc[i] = 0.0;
         }

         // Loop over monomer types (k is a monomer index)
         for (k = 0; k < nMonomer; ++k) {
            double vec = simulator().chiEvecs(j, k)/double(nMonomer);
            // Loop over grid points
            RField<D> const & Wr = w0_[k];
            for (i = 0; i < meshSize; ++i) {
               Wc[i] += vec*Wr[i];
            }
         }
      }
      
      #if 0
      Log::file() << "wc " << wc0_.capacity() << "\n";
      for (i = 0; i < 10; ++i) {
         Log::file() << "wc_1 " << wc0_[0][i] << "\n";
         Log::file() << "wc_2 " << wc0_[1][i] << "\n";
      }
      #endif
   }
   
}
}
#endif 

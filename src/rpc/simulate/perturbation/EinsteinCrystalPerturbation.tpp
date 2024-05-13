#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP

#include "EinsteinCrystalPerturbation.h"
#include <rpc/simulate/Simulator.h>
#include <rpc/System.h>
#include <prdc/cpu/RField.h>
#include <util/containers/DArray.h>
#include <util/global.h>
#include <util/format/Dbl.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(Simulator<D>& simulator)
    : Perturbation<D>(simulator),
      lambda_(0.0),
      dLambda_(0.0),
      f_(0.0),
      df_(0.0),
      interval_(0),
      hamiltonianEC_(0.0),
      hamiltonianBCP_(0.0)
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
      read(in, "outputFileName", outputFileName_);
      read(in,"lambda", lambda_);
      read(in,"rampingRate", dLambda_);
      read(in, "outputInterval", interval_);
      system().fileMaster().openOutputFile(outputFileName_, outputFile_);
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
      double lambda = lambda_ + simulator().iStep()* dLambda_;
      
      // Compute EC Hamiltonian
      double prefactor, w, s;
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      hamiltonianEC_ = 0;
      for (int j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & Wc = simulator().wc(j);
         prefactor = -0.5*double(nMonomer)/simulator().chiEval(j);
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
  
      return lambda* hamiltonianBCP_ + (1.0-lambda) * hamiltonianEC_;
   }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void 
   EinsteinCrystalPerturbation<D>::modifyDc(DArray< RField<D> > & dc)
   {}
   
   
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
   
   /*
   * Compute and output the derivative of the free energy with respect to lambda
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::updateDf(){
      df_ = hamiltonianBCP_ - hamiltonianEC_;
      f_ += df_ * dLambda_;
      double lambda = lambda_ + simulator().iStep() * dLambda_;
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      if (simulator().iStep() % interval_ == 0) {
         outputFile_<< Dbl(lambda);
         outputFile_<< Dbl(f_/nMonomerSystem);
         outputFile_<< "   ";
         outputFile_ << "\n";
      }
   
   }

}
}
#endif 

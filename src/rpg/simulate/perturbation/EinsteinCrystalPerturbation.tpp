#ifndef RPG_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RPG_EINSTEIN_CRYSTAL_PERTURBATION_TPP

#include "EinsteinCrystalPerturbation.h"
#include <rpg/simulate/Simulator.h>
#include <rpg/System.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/Field.h>
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(Simulator<D>& simulator)
    : Perturbation<D>(simulator),
      alpha_(0.0),
      ecHamiltonian_(0.0),
      bcpHamiltonian_(0.0),
      stateEcHamiltonian_(0.0),
      stateBcpHamiltonian_(0.0)
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
      read(in, "lambda", lambda_);
      read(in, "alpha", alpha_);
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
      
      // Allocate memory for reference field
      w0_.allocate(nMonomer);
      wc0_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0_[i].allocate(dimensions);
         wc0_[i].allocate(dimensions);
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
   double 
   EinsteinCrystalPerturbation<D>::hamiltonian(double unperturbedHamiltonian)
   {
      // Compute Einstein crystal Hamiltonian
      double prefactor, s;
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      RField<D> wcs;
      wcs.allocate(meshSize);
      ecHamiltonian_ = 0.0;
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      for (int j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & Wc = simulator().wc(j);
         prefactor = alpha_;
         s = simulator().sc(j);
         assignRealSubtractDouble<<<nBlocks, nThreads>>>
            (wcs.cField(), Wc.cField(), s, meshSize);
         pointWiseSubtract<<<nBlocks, nThreads>>>
            (wcs.cField(), wc0_[j].cField(), meshSize);
         double wSqure = 0;
         wSqure = 
              (double)gpuInnerProduct(wcs.cField(), 
                                      wcs.cField(), meshSize);
         ecHamiltonian_ += prefactor * wSqure;
      }
      
      // Normalize EC hamiltonian to equal a value per monomer
      ecHamiltonian_ /= double(meshSize);
      
      // Compute EC hamiltonian of system
      ecHamiltonian_ *= nMonomerSystem;
      
      // Obtain block copolymer hamiltonian
      bcpHamiltonian_ = unperturbedHamiltonian;
      
      // Compute perturbation to Hamiltonian
      double compH;
      compH = lambda_* bcpHamiltonian_ + (1.0-lambda_) * ecHamiltonian_; 
      return compH - unperturbedHamiltonian;
   }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void 
   EinsteinCrystalPerturbation<D>::incrementDc(DArray< RField<D> > & dc)
   {
      double prefactor, s, b;
      const int meshSize = system().domain().mesh().size();
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();   
      RField<D> DcBCP, DcEC, wRef;
      DcBCP.allocate(meshSize);
      DcEC.allocate(meshSize);
      wRef.allocate(meshSize);
      b = 1.0;
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         RField<D>& Dc = dc[i];
         RField<D> const & Wc = simulator().wc(i);
         s = simulator().sc(i);
         prefactor = -1.0*double(nMonomer)/simulator().chiEval(i)/vMonomer;
         
         // Copy block copolymer derivative
         assignReal<<<nBlocks,nThreads>>>(DcBCP.cField(), 
                                          Dc.cField(), meshSize);
         // Substract Reference field 
         assignReal<<<nBlocks,nThreads>>>(wRef.cField(), 
                                          wc0_[i].cField(), meshSize);
         scaleReal<<<nBlocks,nThreads>>>(wRef.cField(), -1.0, meshSize);
         // Compute EC derivative
         computeDField<<<nBlocks, nThreads>>>
            (DcEC.cField(), Wc.cField(), wRef.cField(), prefactor, b, s, meshSize);
       
         // Compute composite derivative
         scaleReal<<<nBlocks, nThreads>>>(Dc.cField(), lambda_, meshSize);
         pointWiseAddScale<<<nBlocks, nThreads>>>(Dc.cField(), DcEC.cField(), 1.0 - lambda_, meshSize);
      }
   }
   
   /*
   * Compute and return derivative of free energy with respect to lambda.
   */ 
   template <int D>
   double EinsteinCrystalPerturbation<D>::df(){
      return bcpHamiltonian_ - ecHamiltonian_;
   }
   
   /*
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::saveState(){
      stateEcHamiltonian_ = ecHamiltonian_;
      stateBcpHamiltonian_ = bcpHamiltonian_;
   }
   
   /*
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::restoreState(){
      ecHamiltonian_ = stateEcHamiltonian_;
      bcpHamiltonian_ = stateBcpHamiltonian_;
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
      int i, j;
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      DArray<RField<D>> const * Wr = &w0_;
      // Loop over eigenvectors (i is an eigenvector index)
      for (i = 0; i < nMonomer; ++i) {

         // Loop over grid points to zero out field wc0_[i]
         RField<D>& Wc = wc0_[i];
         assignUniformReal<<<nBlocks, nThreads>>>(Wc.cField(), 0, meshSize);
         
         // Loop over monomer types (j is a monomer index)
         for (j = 0; j < nMonomer; ++j) {
            cudaReal vec;
            vec = (cudaReal) simulator().chiEvecs(i, j)/nMonomer;

            // Loop over grid points
            pointWiseAddScale<<<nBlocks, nThreads>>>
               (Wc.cField(), (*Wr)[j].cField(), vec, meshSize);
         }
         
      }
      
      // Debugging output
      std::string filename = "wc";
      system().fieldIo().writeFieldsRGrid(filename, wc0_, 
                                          system().domain().unitCell());
   }

}
}
#endif 

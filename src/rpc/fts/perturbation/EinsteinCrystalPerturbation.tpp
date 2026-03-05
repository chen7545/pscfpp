#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP

#include "EinsteinCrystalPerturbation.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(
                                                  Simulator<D>& simulator)
    : Perturbation<D>(simulator),
      ecHamiltonian_(0.0),
      unperturbedHamiltonian_(0.0),
      stateEcHamiltonian_(0.0),
      stateUnperturbedHamiltonian_(0.0),
      hasEpsilon_(false)
   {  ParamComposite::setClassName("EinsteinCrystalPerturbation"); }

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
      ParamComposite::read(in, "lambda", lambda_);

      // Allocate and initialize epsilon_ array
      const int nMonomer = system().mixture().nMonomer();
      epsilon_.allocate(nMonomer - 1);
      for (int i = 0; i < nMonomer - 1 ; ++i) {
         epsilon_[i] = 0.0;
      }

      // Optionally read the parameters used in Einstein crystal integration
      int nc = nMonomer - 1;
      hasEpsilon_
         = readOptionalDArray(in, "epsilon", epsilon_, nc).isActive();

      ParamComposite::read(in, "fieldFileName", fieldFileName_);
   }

   /*
   * Setup before simulation.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::setup()
   {
      const int nMonomer = system().mixture().nMonomer();
      const IntVec<D> meshDimensions
                      = system().domain().mesh().dimensions();

      // Allocate memory for reference field
      w0_.allocate(nMonomer);
      wc0_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0_[i].allocate(meshDimensions);
         wc0_[i].allocate(meshDimensions);
      }
      dw_.allocate(meshDimensions);

      /*
      * If the user did not input epsilon_ values, set values to -1.0
      * times the nontrivial eigenvalues of the projected chi matrix.
      */
      if (!hasEpsilon_){
         for (int i = 0; i < nMonomer - 1 ; ++i) {
            epsilon_[i] = -1.0 * simulator().chiEval(i);
         }
      }

      // Check that all epsilon values are positive
      for (int i = 0; i < nMonomer - 1 ; ++i) {
         UTIL_CHECK(epsilon_[i] > 0.0);
      }

      // Read in reference field from a file
      UnitCell<D> tempUnitCell;
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGrid(fieldFileName_, w0_, tempUnitCell);

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
      // Set unperturbedHamiltonian_ member variable
      unperturbedHamiltonian_ = unperturbedHamiltonian;

      // Constants
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      Domain<D> const & domain = system().domain();
      const int meshSize = domain.mesh().size();
      const double vSystem  = domain.unitCell().volume();
      const double nMonomerSystem = vSystem / vMonomer;

      // Compute Einstein crystal Hamiltonian
      double prefactor;
      ecHamiltonian_ = 0.0;
      for (int j = 0; j < nMonomer - 1; ++j) {
         VecOp::subVV(dw_, simulator().wc(j), wc0_[j]);
         prefactor = double(nMonomer)/(2.0 * epsilon_[j]);
         ecHamiltonian_ += prefactor * Reduce::sumSq(dw_);
      }
      ecHamiltonian_ /= double(meshSize);
      ecHamiltonian_ *= nMonomerSystem;

      return (1.0 - lambda_)*(ecHamiltonian_ - unperturbedHamiltonian_);
   }

   /*
   * Modify functional derivatives.
   */
   template <int D>
   void
   EinsteinCrystalPerturbation<D>::incrementDc(DArray< RField<D> > & dc)
   {
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      double prefactor;

      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         prefactor = double(nMonomer) / (epsilon_[i] * vMonomer);

         // dw_ = prefactor * (wc - wc0_)
         VecOp::subVV(dw_, simulator().wc(i), wc0_[i]);
         VecOp::mulEqS(dw_, prefactor);

         // dc = dc * lambda_ + dw_ * (1 - lambda_)
         VecOp::mulEqS(dc[i], lambda_);
         VecOp::addEqVc(dc[i], dw_, 1.0 - lambda_);
      }
   }

   /*
   * Compute and return derivative of free energy with respect to lambda.
   */
   template <int D>
   double EinsteinCrystalPerturbation<D>::df()
   {  return unperturbedHamiltonian_ - ecHamiltonian_; }

   /*
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::saveState()
   {
      stateEcHamiltonian_ = ecHamiltonian_;
      stateUnperturbedHamiltonian_ = unperturbedHamiltonian_;
   }

   /*
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::restoreState()
   {
      ecHamiltonian_ = stateEcHamiltonian_;
      unperturbedHamiltonian_ = stateUnperturbedHamiltonian_;
   }

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs of the projected chi matrix as a basis.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::computeWcReference()
   {
      const int nMonomer = system().mixture().nMonomer();
      int i, j;

      // Loop over eigenvectors
      for (i = 0; i < nMonomer; ++i) {

         // Initialize to zero
         VecOp::eqS(wc0_[i], 0.0);

         // Loop over monomer types
         for (j = 0; j < nMonomer; ++j) {
            double vec = simulator().chiEvecs(i, j)/double(nMonomer);
            VecOp::addEqVc(wc0_[i], w0_[j], vec);
         }
      }
   }

}
}
#endif

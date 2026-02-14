#ifndef RPC_FOURTH_ORDER_PARAMETER_TPP
#define RPC_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/ReduceCx.h>




#include <util/containers/DArray.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D>::FourthOrderParameter(Simulator<D>& simulator,
                                                 System<D>& system)
    : AverageAnalyzer<D>(simulator, system),
      kSize_(1),
      isInitialized_(false)
   {  ParamComposite::setClassName("FourthOrderParameter"); }

   /*
   * Destructor.
   */
   template <int D>
   FourthOrderParameter<D>::~FourthOrderParameter()
   {}

   /*
   * Setup before the main loop.
   */
   template <int D>
   void FourthOrderParameter<D>::setup()
   {
      // Precondition: The system must have exactly two monomer types
      UTIL_CHECK(system().mixture().nMonomer() == 2);

      AverageAnalyzer<D>::setup();

      // Compute k-space mesh kMeshDimensions_ and kSize_
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wK_.allocate(dimensions);
         prefactor_.allocate(kMeshDimensions_);
         psi_.allocate(kMeshDimensions_);
         VecOp::eqS(prefactor_, 0.0);
      }

      computePrefactor();
      isInitialized_ = true;
   }

   /*
   * Compute and return the order parameter.
   */
   template <int D>
   double FourthOrderParameter<D>::compute()
   {
      UTIL_CHECK(isInitialized_);
      UTIL_CHECK(wK_.capacity() == kSize_);
      UTIL_CHECK(prefactor_.capacity() == kSize_);
      UTIL_CHECK(psi_.capacity() == kSize_);
      UTIL_CHECK(system().w().hasData());

      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      // Fourier transform W_(r) to obtain wK_ = W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);

      // Evaluate fourth powers, scaled by prefactors
      VecOp::sqSqAbsV(psi_, wK_);
      VecOp::mulEqV(psi_, prefactor_);

      // Summation
      double orderParameter = Reduce::sum(psi_, 1, kSize_);
      orderParameter = std::pow(orderParameter, 0.25);

      return orderParameter;
   }

   /*
   * Output a sampled or block average value.
   */
   template <int D>
   void FourthOrderParameter<D>::outputValue(int step, double value)
   {
      int nSamplePerOutput = AverageAnalyzer<D>::nSamplePerOutput();
      if (simulator().hasRamp() && nSamplePerOutput == 1) {
         std::ofstream& file = AverageAnalyzer<D>::outputFile_;
         UTIL_CHECK(file.is_open());
         double chi = system().interaction().chi(0,1);
         file << Int(step);
         file << Dbl(chi);
         file << Dbl(value);
         file << "\n";
      } else {
         AverageAnalyzer<D>::outputValue(step, value);
      }
   }

   /*
   * Compute prefactors for all wavevectors.
   */
   template <int D>
   void FourthOrderParameter<D>::computePrefactor(Array<double>& prefactor)
   {
      IntVec<D> G;
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      DArray<IntVec<D>> GminList;
      GminList.allocate(kSize_);
      MeshIterator<D> itr(kMeshDimensions_);
      MeshIterator<D> searchItr(kMeshDimensions_);

      // Calculate GminList
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         GminList[itr.rank()] = Gmin;
      }

      // Compute prefactor for each wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If prefactor of current wavevector has not been assigned
         if (prefactor[itr.rank()] == 0){
            Gmin = GminList[itr.rank()];

            // Compute inverse of wavevector
            nGmin.negate(Gmin);

            // Search for inverse of wavevector
            searchItr = itr;
            for (; !searchItr.atEnd(); ++searchItr){
               if (nGmin == GminList[searchItr.rank()]){
                  prefactor[itr.rank()] = 1.0/2.0;
                  prefactor[searchItr.rank()] = 1.0/2.0;
                  inverseFound = true;
               }
            }

            if (inverseFound == false){
               prefactor[itr.rank()]  = 1.0;
            }

         }

      }
   }

   /*
   * Initialize prefactor_ member variable.
   */
   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
   {  computePrefactor(prefactor_); }

}
}
#endif

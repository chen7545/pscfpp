#ifndef RPG_FOURTH_ORDER_PARAMETER_TPP
#define RPG_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>

#include <prdc/cuda/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/cudaTypes.h>
#include <pscf/cpu/VecOp.h>

#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

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
   * FourthOrderParameter setup
   */
   template <int D>
   void FourthOrderParameter<D>::setup()
   {
      // Local copies of data
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Precondition: Require that the system has two monomer types
      UTIL_CHECK(nMonomer == 2);

      AverageAnalyzer<D>::setup();

      // Compute DFT k-space mesh kMeshDimensions_ and kSize_
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

      // Convert W_(r) to fourier mode W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);

      // psi_[i] = |wK_[i]|^4 * prefactor[i]
      VecOp::sqSqAbsV(psi_, wK_);
      VecOp::mulEqV(psi_, prefactor_);

      // Summation
      double orderParameter = Reduce::sum(psi_, 1, kSize_);
      orderParameter = std::pow(orderParameter, 0.25);

      return orderParameter;
   }

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

      // Compute weight factor for each G wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If the weight factor of the current wavevector has not been assigned
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

   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
   {
      HostDArray<cudaReal> prefactor_h(kSize_);
      VecOp::eqS(prefactor_h, 0.0);

      computePrefactor(prefactor_h);

      // Copy the weight factor from cpu(host) to gpu(device)
      prefactor_ = prefactor_h;
   }


}
}
#endif

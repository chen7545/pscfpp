#ifndef RPC_FOURTH_ORDER_PARAMETER_TPP
#define RPC_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/ReduceCx.h>

#include <util/containers/DArray.h>
#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>
#include <complex>
#include <vector>
#include <numeric>
#include <cmath>
#include <set>

namespace Pscf {
namespace Rpc {

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
      AverageAnalyzer<D>::setup();

      // Precondition: Require that the system has two monomer types
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wK_.allocate(dimensions);
         prefactor_.allocate(kSize_);
         psi_.allocate(kSize_);
         VecOp::eqS(prefactor_, 0.0);
      }

      isInitialized_ = true;

      computePrefactor();
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

      // Transform W_(r) to fourier representation W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);

      // Evaluate sum of fourth powers, scaled by array of prefactors
      VecOp::sqSqAbsV(psi_, wK_);
      VecOp::mulEqV(psi_, prefactor_);
      double orderParameter = Reduce::sum(psi_, 1, kSize_);

      orderParameter = std::pow(orderParameter, 0.25);
      return orderParameter;

      #if 0
      // Debugging output
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> G;
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      double kSq;
      std::vector<double> k(kSize_);

      // Calculate GminList
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         kSq = unitCell.ksq(Gmin);
         k[itr.rank()] = kSq;
      }

      auto maxIt = std::max_element(psi.begin(), psi.end());

      // Calculate the index of the maximum element
      size_t maxIndex = std::distance(psi.begin(), maxIt);
      double kmax = k[maxIndex];

      Log::file() << std::endl;
      for (itr.begin(); !itr.atEnd(); ++itr){
         if (k[itr.rank()] == kmax){
            G = itr.position();
            Gmin = shiftToMinimum(G, meshDimensions, unitCell);
            Log::file() << "ksq: " << k[itr.rank()] << std::endl;
            Log::file() << " G: " << G<< std::endl;
            Log::file() << " Gmin: " << Gmin<< std::endl;
            Log::file() << " prefactor: " <<  prefactor_[itr.rank()]<< std::endl;
            Log::file() << " psi: " <<  psi[itr.rank()]<< std::endl;
         }

      }
      #endif

   }

   template <int D>
   void FourthOrderParameter<D>::outputValue(int step, double value)
   {
      int nSamplePerOutput = AverageAnalyzer<D>::nSamplePerOutput();
      if (simulator().hasRamp() && nSamplePerOutput == 1) {
         std::ofstream& file = AverageAnalyzer<D>::outputFile_;
         UTIL_CHECK(file.is_open());
         double chi= system().interaction().chi(0,1);
         file << Int(step);
         file << Dbl(chi);
         file << Dbl(value);
         file << "\n";
      } else {
         AverageAnalyzer<D>::outputValue(step, value);
      }
   }

   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
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

      // Compute prefactor for each G wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If prefactor of current wavevector has not been assigned
         if (prefactor_[itr.rank()] == 0){
            Gmin = GminList[itr.rank()];

            // Compute inverse of wavevector
            nGmin.negate(Gmin);

            // Search for inverse of wavevector
            searchItr = itr;
            for (; !searchItr.atEnd(); ++searchItr){
               if (nGmin == GminList[searchItr.rank()]){
                  prefactor_[itr.rank()] = 1.0/2.0;
                  prefactor_[searchItr.rank()] = 1.0/2.0;
                  inverseFound = true;
               }
            }

            if (inverseFound == false){
               prefactor_[itr.rank()]  = 1.0;
            }

         }

      }
   }

}
}
#endif

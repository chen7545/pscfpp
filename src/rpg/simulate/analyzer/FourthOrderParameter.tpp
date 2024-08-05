#ifndef RPG_FOURTH_ORDER_PARAMETER_TPP
#define RPG_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpg/simulate/Simulator.h>
#include <rpg/System.h>

#include <prdc/cuda/RField.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <fftw3.h>

#include <iostream>
#include <complex>
#include <vector>
#include <numeric>
#include <cmath>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D>::FourthOrderParameter(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      kSize_(1),
      hasAverage_(true),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("FourthOrderParameter"); }


   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void FourthOrderParameter<D>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOptional(in, "hasAverage", hasAverage_);
      readOutputFileName(in);
      readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
      
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "    chi       " << "FourthOrderParameter" << "\n";
   }
   
   /*
   * FourthOrderParameter setup
   */
   template <int D>
   void FourthOrderParameter<D>::setup() 
   {
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      if (nMonomer != 2) {
         UTIL_THROW("The FourthOrderParameter Analyzer is designed specifically for diblock copolymer system. Please verify the number of monomer types in your system.");
      }
      
      //Allocate variables
      IntVec<D> const & dimensions = system().mesh().dimensions();
      if (!isInitialized_){
         wK_.allocate(dimensions);
      }
      
      // Compute Fourier space dimension
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
            kSize_ *= dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
            kSize_ *= (dimensions[i]/2 + 1);
         }
      }
      
      isInitialized_ = true;
      
      // Clear accumulators
      if (hasAverage_){
         accumulator_.clear();
      }
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }

   /* 
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void FourthOrderParameter<D>::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;
      computeFourthOrderParameter();
      
      if (hasAverage_){
         accumulator_.sample(FourthOrderParameter_);
      }
      
      double chi =  system().interaction().chi(0,1);
      UTIL_CHECK(outputFile_.is_open());
      outputFile_ << Dbl(chi);
      outputFile_ << Dbl(FourthOrderParameter_);
      outputFile_<< "\n";
   }
   
   template <int D>
   void FourthOrderParameter<D>::computeFourthOrderParameter()
   {
      RField<D> psi;
      psi.allocate(kSize_);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(kSize_, nBlocks, nThreads);
      
      // Conver W_(r) to fourier mode W_(k)
      system().fft().forwardTransform(simulator().wc(0), wK_);
      
      // W_(k)^2
      squaredMagnitudeComplex<<<nBlocks,nThreads>>>
            (wK_.cField(), psi.cField(), kSize_);
      
      // W_(k)^4
      inPlacePointwiseMul<<<nBlocks,nThreads>>>
            (psi.cField(), psi.cField(), kSize_);
            
      // Get sum over all wavevectors
      FourthOrderParameter_ = (double)gpuSum(psi.cField(), kSize_);
      FourthOrderParameter_ = std::pow(FourthOrderParameter_, 0.25);
   }
   
   /*
   * Output final results to output file.
   */
   template <int D>  
   void FourthOrderParameter<D>::output() 
   {
      if (hasAverage_){
         Log::file() << std::endl;
         Log::file() << "At chi = " << system().interaction().chi(0,1) << "\n";
         Log::file() << "Time average of the FourthOrderParameter is: "
                     << Dbl(accumulator_.average())
                     << " +- " << Dbl(accumulator_.blockingError(), 9, 2) 
                     << "\n";
      } 
      
   }

}
}
#endif 

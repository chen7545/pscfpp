#ifndef RPC_CHI_NORMALIZATION_TPP
#define RPC_CHI_NORMALIZATION_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiNormalization.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/System.h>

#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>

#include <prdc/cpu/RField.h>
#include <prdc/crystal/UnitCell.h>

#include <util/format/Int.h>
#include <util/misc/Log.h>
#include <cmath> 

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ChiNormalization<D>::ChiNormalization(Simulator<D>& simulator, System<D>& system)
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      isInitialized_(false),
      nBond_(100)
   {  setClassName("ChiNormalization"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ChiNormalization<D>::readParameters(std::istream& in) 
   { 
      readOutputFileName(in); 
   }
   
   /*
   * ChiNormalization setup
   */
   template <int D>
   void ChiNormalization<D>::setup() 
   {
      // Only work for conformationally symmetric polymer
      double kuhnA = system().mixture().polymer(0).block(0).kuhn();
      double kuhnB = system().mixture().polymer(0).block(1).kuhn();
      UTIL_CHECK(std::abs(kuhnA/kuhnB -1) < 1e-5);
      
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      
      #if 0
      UTIL_CHECK(unitCell.lattice() == UnitCell<3>::Cubic || 
                 unitCell.lattice() == UnitCell<3>::Tetragonal || 
                 unitCell.lattice() == UnitCell<3>::Orthorhombic)
      #endif
      double nA = system().mixture().polymer(0).block(0).length();
      double nB = system().mixture().polymer(0).block(1).length();
      int N = nA + nB;
      for (int i = 0; i < D; i++){
         // The box length with unit R0 (a*sqrt(N))
         dx_[i] = unitCell.rBasis(i)[i] / dimensions[i] * sqrt(N);
      }
      
   }
   
   /**
   * Return and compute chi linear normalization term.
   */
   template <int D>
   double ChiNormalization<D>::computeLinearNormalization()
   {
      const double vMonomer = system().mixture().vMonomer();
      
      double vCell = 1.0;
      for (int i = 0; i < D; i++){
         vCell *= dx_[i];
      }
      double sumP = 0;
      
      for (int i = 1; i < nBond_; i++){
         sumP += computeP(i);
      }
      
      // Add contribution i > nBond (For large n, erf(f(n)) close to 1)
      sumP += (2.0 / std::sqrt(0.5 + nBond_))* std::pow(3.0 / (2.0 * M_PI), 1.5) * vCell;

      double zInf = 1.0 - (1.0 + 2.0 * sumP) / vCell*vMonomer;
      return zInf;
   }
   
   /**
   * Return and compute the probability function.
   */ 
   template <int D>
   double ChiNormalization<D>::computeP(int i)
   {
      double p = 1.0;
      for (int d = 0; d < D; d++){
         p *= dx_[d] * sqrt(3.0 / (2.0 * M_PI * i)) * erf(M_PI / dx_[d] * sqrt(i/6.0)); 
      }
      return p;
   }

   /*
   * Output final results to output file.
   */
   template <int D>  
   void ChiNormalization<D>::output() 
   {
      zInf_ = computeLinearNormalization();
      
      Log::file() << "Chi linear normalization (zInf) is : "<< zInf_;

      // Output linear normalization to one file
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "\t" << "Chi linear normalization (zInf) is : ";
      outputFile_ << zInf_;
      outputFile_.close();
   }

}
}
#endif

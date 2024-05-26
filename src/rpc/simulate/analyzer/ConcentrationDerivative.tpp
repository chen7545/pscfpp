#ifndef RPC_CONCENTRATION_DERIVATIVE_TPP
#define RPC_CONCENTRATION_DERIVATIVE_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"
#include "Analyzer.h"
#include <rpc/simulate/Simulator.h>
#include <rpc/System.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <sstream>

namespace Pscf {
namespace Rpc 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ConcentrationDerivative<D>::ConcentrationDerivative(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      nSamplePerBlock_(1),
      isInitialized_(false),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {  setClassName("ConcentrationDerivative"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ConcentrationDerivative<D>::readParameters(std::istream& in) 
   {
      Analyzer<D>::readParameters(in);
      
      // Read parameters for accumulator
      read(in,"nSamplePerBlock", nSamplePerBlock_);
      
      isInitialized_ = true;
   }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ConcentrationDerivative<D>::setup() 
   {  
      std::string filename;
      filename_  = outputFileName();
      system().fileMaster().openOutputFile(filename_ , outputFile_);
      
      accumulator_.setNSamplePerBlock(nSamplePerBlock_);
      accumulator_.clear();
   }

   /*
   * Periodically write a frame to file
   */
   template <int D>
   void ConcentrationDerivative<D>::sample(long iStep) 
   {  
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2); 
      
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();
      
      if (!isAtInterval(iStep)) return;
      
      // Obteain Hamiltonian per monomer
      double h = simulator().hamiltonian()/nMonomerSystem;
      
      // Calculate derivative with respect to concentration
      double dfdc = h * vMonomer;
      
      // With N term
      double Hh = meshSize/2/nMonomerSystem;
      //dfdc -= Hh;
      accumulator_.sample(dfdc);
      
      if (nSamplePerBlock_ > 0) { 
         if (accumulator_.isBlockComplete()) {
            UTIL_CHECK(outputFile_.is_open());
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep);
            double block = accumulator_.blockAverage();
            outputFile_ << Dbl(block);
         }
            outputFile_ << "\n";
      }
   }
  
   /*
   * Close output file at end of simulation.
   */
   template <int D>
   void ConcentrationDerivative<D>::output() 
   {   // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
      
      // Write average (*.ave) and error analysis (*.aer) file
      outputAccumulator(outputFileName());
  
   }
   
   
   /*
   * Output results to file after simulation is completed.
   */
   template <int D> 
   void ConcentrationDerivative<D>::outputAccumulator(std::string outputFileName)
   {
   
      // Close data (*.dat) file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      // Write average (*.ave) file
      std::string fileName;
      fileName = outputFileName;
      fileName += ".ave";
      system().fileMaster().openOutputFile(fileName, outputFile_);
      double ave, err;
      ave = accumulator_.average();
      err = accumulator_.blockingError();
      outputFile_ << " " << std::left << std::setw(2) 
                      << "dH/dchi" << "   ";
      outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
                      
      outputFile_.close();

      // Write error analysis (*.aer) file
      fileName = outputFileName;
      fileName += ".aer";
      system().fileMaster().openOutputFile(fileName, outputFile_);
      std::string line; 
      line = 
      "---------------------------------------------------------------------";
      outputFile_ << line << std::endl;
      accumulator_.output(outputFile_);
      outputFile_ << std::endl;
      outputFile_.close();
   
   }

}
}
#endif

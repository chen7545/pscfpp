#ifndef RPC_CHI_DERIVATIVE_TPP
#define RPC_CHI_DERIVATIVE_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"
#include "Analyzer.h"
#include <rpc/simulate/Simulator.h>
#include <rpc/System.h>
#include <util/misc/FileMaster.h>
//#include <util/archives/Serializable_includes.h>
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
   ChiDerivative<D>::ChiDerivative(Simulator<D>& simulator, 
                                   System<D>& system) 
    : Analyzer<D>(),
      nSamplePerBlock_(1),
      isInitialized_(false),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {  setClassName("ChiDerivative"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ChiDerivative<D>::readParameters(std::istream& in) 
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
   void ChiDerivative<D>::setup() 
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
   void ChiDerivative<D>::sample(long iStep) 
   {  
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2); 
      
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();
      double chi= system().interaction().chi(0,1);
      
      if (!isAtInterval(iStep)) return;
      
      // Obteain fieldHamiltonian per monomer
      double hField = simulator().fieldHamiltonian()/nMonomerSystem;

      // Compute current derivative of free energy with repect to chi_bare
      double dfdchi = -(hField - 0.5*simulator().sc(nMonomer - 1))/chi + 1.0/4.0;
      
      // With N term
      //dfdchi += meshSize/(2* nMonomerSystem* chi);
      accumulator_.sample(dfdchi*nMonomerSystem);
      
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
   void ChiDerivative<D>::output() 
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
   void ChiDerivative<D>::outputAccumulator(std::string outputFileName)
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

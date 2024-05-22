#ifndef RPC_PERTURBATION_DERIVATIVE_TPP
#define RPC_PERTURBATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"
#include <rpc/System.h>
#include <rpc/simulate/Simulator.h>
#include <rpc/simulate/perturbation/Perturbation.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   PerturbationDerivative<D>::PerturbationDerivative(Simulator<D>& simulator, 
                                               System<D>& system)
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      outputInterval_(0),
      nSamplePerBlock_(0),
      f_(0.0),
      isInitialized_(false),
      hasAccumulator_(false)
   {  setClassName("PerturbationDerivative"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void PerturbationDerivative<D>::readParameters(std::istream& in) 
   {
      // Read sample interval
      readInterval(in);
   
      // Read output interval and file name
      readOptional(in, "outputInterval", outputInterval_);
      readOutputFileName(in);
      
      // Read parameters for autocorrelation
      readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
      
      // Open outputFile
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
   }

   /*
   * PerturbationDerivative setup
   */
   template <int D>
   void PerturbationDerivative<D>::setup() 
   {
      // Set all parameters and allocate to initialize state for autocorrelation 
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      if (simulator().perturbation().mode() == 0){
         UTIL_CHECK(nSamplePerBlock_ > 0);
      }
      if (simulator().perturbation().mode() == 1){
         UTIL_CHECK(outputInterval_ > 0);
      }
      accumulator_.setNSamplePerBlock(nSamplePerBlock_);
      accumulator_.clear();
      hasAccumulator_ = true;
      isInitialized_ = true;
   }
   
   /*
   * Sample hamiltonian and add to accumulator.
   */
   template <int D>
   void PerturbationDerivative<D>::sample(long iStep) 
   {
      UTIL_CHECK(simulator().hasPerturbation());
      if (!isAtInterval(iStep)) return;
      updateAccumulator(iStep);
   }
   
   template <int D>
   void PerturbationDerivative<D>::output()
   {
      
      UTIL_CHECK(hasAccumulator_);

      // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
      
      // Write average (*.ave) and error analysis (*.aer) file
      if (simulator().perturbation().mode() == 0){
         outputAccumulator(outputFileName());
      }
      
   }

   /*
   * Update accumulator for current value.
   */
   template <int D>   
   void PerturbationDerivative<D>::updateAccumulator(long iStep) 
   {
      UTIL_CHECK(hasAccumulator_);
      
      double df = simulator().perturbation().df();
      
      // Static parameter
      if (simulator().perturbation().mode() == 0){
         
         // Update accumulators
         accumulator_.sample(df);
         
         // Output block averages
         if (accumulator_.isBlockComplete()) {
            UTIL_CHECK(outputFile_.is_open());
            outputFile_ << Dbl(simulator().perturbation().lambda());
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep);
            double block = accumulator_.blockAverage();
            outputFile_ << Dbl(block);
            outputFile_ << "\n";
         }
         
      } else if (simulator().perturbation().mode() == 1){
         
         //Output countinuous thermodynamic integration
         if(iStep % outputInterval_ == 0){
            UTIL_CHECK(outputFile_.is_open());
            outputFile_ << Dbl(simulator().perturbation().lambda());
            outputFile_ << Dbl(f_);
            outputFile_ << "\n";
         }
         
         // Update free energy
         f_ += df * interval();
      }
   } 
   
   /*
   * Output results to file after simulation is completed.
   */
   template <int D> 
   void PerturbationDerivative<D>::outputAccumulator(std::string outputFileName)
   {
      UTIL_CHECK(hasAccumulator_);

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
                      << "df/dX" << "   ";
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

#ifndef RPG_AM_ERROR_ANALYZER_TPP
#define RPG_AM_ERROR_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmErrorAnalyzer.h"

#include <rpg/System.h>
#include <rpg/simulate/Simulator.h>
#include <rpg/simulate/compressor/Compressor.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Rpg 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   AmErrorAnalyzer<D>::AmErrorAnalyzer(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      nSamplePerBlock_(0),
      nValue_(100),
      hasAccumulators_(false),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {setClassName("AmErrorAnalyzer");}

   /*
   * Destructor.
   */
   template <int D>
   AmErrorAnalyzer<D>::~AmErrorAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void AmErrorAnalyzer<D>::readParameters(std::istream& in) 
   {
      Analyzer<D>::readParameters(in);
      readOptional(in, "AmMaxItr", nValue_);
      read(in,"nSamplePerBlock", nSamplePerBlock_);
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      projectionRatioAccumulators_.allocate(nValue_);
      mixingRatioAccumulators_.allocate(nValue_);
      predictRatioAccumulators_.allocate(nValue_);
      for (int i = 0; i < nValue_; ++i) {
         projectionRatioAccumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
         projectionRatioAccumulators_[i].clear();
         mixingRatioAccumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
         mixingRatioAccumulators_[i].clear();
         predictRatioAccumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
         predictRatioAccumulators_[i].clear();
      }
   }
   
   /*
   * Setup before system.
   */ 
   template <int D>
   void AmErrorAnalyzer<D>::setup()
   {
      
      for (int i = 0; i < nValue_; ++i) {
         projectionRatioAccumulators_[i].clear();
         mixingRatioAccumulators_[i].clear();
         predictRatioAccumulators_[i].clear();
      }
      
      hasAccumulators_ = true;
   }

   /*
   * Compute and sample current values.
   */
   template <int D>
   void AmErrorAnalyzer<D>::sample(long iStep) 
   {
      UTIL_CHECK(hasAccumulators_);
      if (!isAtInterval(iStep)) return;
      
      std::vector<double> stepOneRatio = simulator().compressor().stepOneRatioVector();
      std::vector<double> predictRatio = simulator().compressor().predictRatioVector();
      std::vector<double> stepTwoRatio = simulator().compressor().stepTwoRatioVector();
      // Update accumulators.
      for (int i = 0; i < stepOneRatio.size(); ++i) {
         double data = stepOneRatio[i];
         projectionRatioAccumulators_[i].sample(data);
      }
      
      for (int i = 0; i < predictRatio.size(); ++i) {
         double data = predictRatio[i];
         predictRatioAccumulators_[i].sample(data);
      }
      
      for (int i = 0; i < stepTwoRatio.size(); ++i) {
         double data = stepTwoRatio[i];
         mixingRatioAccumulators_[i].sample(data);
      }

   }

   /*
   * Output results after a system is completed.
   */
   template <int D>
   void AmErrorAnalyzer<D>::output()
   {   
      // Write average (*.ave) and error analysis (*.aer) files
      outputAccumulators(outputFileName());
      
   }
   
   /*
   * Output results to file after simulation is completed.
   */
   template <int D> 
   void AmErrorAnalyzer<D>::outputAccumulators(std::string outputFileName)
   {
      UTIL_CHECK(hasAccumulators_);

      // Close data (*.dat) file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      // Compute maximum length of name field
      int nameWidth = 0;
      int length;
      for (int i = 0; i < nValue_; ++i) {
         length = 3;
         if (length > nameWidth) {
            nameWidth = length;
         }
      }
      nameWidth += 2;

      // Write average (*.ave) file
      std::string fileName;
      fileName = outputFileName;
      fileName += ".ave";
      system().fileMaster().openOutputFile(fileName, outputFile_);
      double aveOne, errOne, aveTwo, errTwo, avePred, errPred;
      outputFile_ << " " << std::left << std::setw(nameWidth) 
                      << "#itr" << "    " << "stepOneRatio" << "  " 
                      << " +- " << "error" << "          " 
                      << "stepTwoRatio" << "  " 
                      << " +- " << "error" << "          " 
                      << "stepOnePredRatio" << "  " 
                      << " +- " << "error" << std::endl;
      for (int i = 0; i < nValue_; ++i) {
         aveOne = projectionRatioAccumulators_[i].average();
         errOne = projectionRatioAccumulators_[i].blockingError();
         avePred = predictRatioAccumulators_[i].average();
         errPred = predictRatioAccumulators_[i].blockingError();
         aveTwo = mixingRatioAccumulators_[i].average();
         errTwo = mixingRatioAccumulators_[i].blockingError();
         outputFile_ << " " << std::left << std::setw(nameWidth) 
                      << i + 1 << "   ";
         outputFile_ << Dbl(aveOne) << " +- " << Dbl(errOne, 6, 2)<< "    ";
         outputFile_ << Dbl(aveTwo) << " +- " << Dbl(errTwo, 6, 2) << "    ";
         outputFile_ << Dbl(avePred) << " +- " << Dbl(errPred, 6, 2)<< "\n";
      }
      outputFile_.close();

      // Write error analysis (*.aer) file
      fileName = outputFileName;
      fileName += ".aer";
      system().fileMaster().openOutputFile(fileName, outputFile_);
      std::string line; 
      line = 
      "---------------------------------------------------------------------";
      for (int i = 0; i < nValue_; ++i) {
         outputFile_ << line << std::endl;
         outputFile_ << "Itr "<< i + 1 << " :" << std::endl;
         projectionRatioAccumulators_[i].output(outputFile_);
         outputFile_ << "   ";
         mixingRatioAccumulators_[i].output(outputFile_);
         outputFile_ << std::endl;
         predictRatioAccumulators_[i].output(outputFile_);
         outputFile_ << "   ";
      }
      outputFile_.close();
      
   }
   
   
   
   
}
}
#endif

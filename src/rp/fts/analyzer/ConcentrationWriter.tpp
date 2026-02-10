#ifndef RP_CONCENTRATION_WRITER_TPP
#define RP_CONCENTRATION_WRITER_TPP
/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationWriter.h"
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   ConcentrationWriter<D,T>::ConcentrationWriter(
                                   typename T::Simulator& simulator,
                                   typename T::System& system)
    : AnalyzerT(simulator, system),
      nSample_(0),
      isInitialized_(false)
   {  ParamComposite::setClassName("ConcentrationWriter"); }

   /*
   * Read interval and outputFileName.
   */
   template <int D, class T>
   void ConcentrationWriter<D,T>::readParameters(std::istream& in)
   {
      AnalyzerT::readParameters(in);
      isInitialized_ = true;
   }

   /*
   * Initialize before main simulation loop.
   */
   template <int D, class T>
   void ConcentrationWriter<D,T>::setup()
   {
      nSample_ = 0;
      std::string filename = AnalyzerT::outputFileName();
      system().fileMaster().openOutputFile(filename, outputFile_);
      writeHeader(outputFile_);
   }

   template <int D, class T>
   void ConcentrationWriter<D,T>::writeFrame(std::ofstream& out, long iStep)
   {
      UTIL_CHECK(system().w().hasData());
      if (!system().c().hasData()){
         system().compute();
      }
      out << "i = " << iStep << "\n";
      bool writeHeader = false;
      bool isSymmetric = false;
      typename T::Domain const & domain = system().domain();
      typename T::FieldIo const & fieldIo = domain.fieldIo();
      fieldIo.writeFieldsRGrid(out, system().c().rgrid(),
                               domain.unitCell(),
                               writeHeader, isSymmetric);
      out << "\n";
   }

   template <int D, class T>
   void ConcentrationWriter<D,T>::writeHeader(std::ofstream& out)
   {
      int nMonomer = system().mixture().nMonomer();
      bool isSymmetric = false;
      typename T::Domain const & domain = system().domain();
      typename T::FieldIo const & fieldIo = domain.fieldIo();
      fieldIo.writeFieldHeader(out, nMonomer,
                               domain.unitCell(), isSymmetric);
      out << "\n";
   }


   /*
   * Periodically write a frame to file
   */
   template <int D, class T>
   void ConcentrationWriter<D,T>::sample(long iStep)
   {
      if (AnalyzerT::isAtInterval(iStep))  {
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
   }

   /*
   * Close output file at end of simulation.
   */
   template <int D, class T>
   void ConcentrationWriter<D,T>::output()
   {  outputFile_.close(); }

}
}
#endif

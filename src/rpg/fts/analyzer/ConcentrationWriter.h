#ifndef RPG_CONCENTRATION_WRITER_H
#define RPG_CONCENTRATION_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"        // Base class template

namespace Pscf {
namespace Rpg {

   // Forward declarations
   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Periodically write snapshots to a trajectory file.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationWriter : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      ConcentrationWriter(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~ConcentrationWriter()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Clear nSample counter.
      */
      virtual void setup();

      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep step index
      */
      virtual void sample(long iStep);

      /**
      * Close trajectory file after run.
      */
      virtual void output();

      using ParamComposite::setClassName;
      using ParamComposite::read;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::system;

   protected:

      // Output file stream
      std::ofstream outputFile_;

      // Output filename
      std::string filename_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;

      /// Has readParam been called?
      long isInitialized_;

   protected:

      /**
      * Write data that should appear once, at beginning of the file.
      *
      * Called by sample on first invocation. Default implementation is empty.
      *
      * \param out output file stream
      */
      void writeHeader(std::ofstream& out);

      /**
      * Write data that should appear in every frame.
      *
      * \param out output file stream
      * \param iStep MC time step index
      */
      void writeFrame(std::ofstream& out, long iStep);

   };

   // Explicit instantiation declarations
   extern template class ConcentrationWriter<1>;
   extern template class ConcentrationWriter<2>;
   extern template class ConcentrationWriter<3>;

}
}
#endif

#ifndef RP_TRAJECTORY_WRITER_H
#define RP_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Periodically write field frames (snapshots) to a trajectory file.
   *
   * \see \ref rp_TrajectoryWriter_page "Manual Page"
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class TrajectoryWriter : public T::Analyzer
   {

   public:

      /**
      * Constructor.
      */
      TrajectoryWriter(typename T::Simulator& simulator, 
                       typename T::System& system);

      /**
      * Destructor.
      */
      virtual ~TrajectoryWriter()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      void readParameters(std::istream& in) override;

      /**
      * Setup before main simulation loop. 
      */
      void setup() override;

      /**
      * Write a frame/snapshot to the trajectory file.
      *
      * \param iStep  step index
      */
      void sample(long iStep) override;

      /**
      * Close trajectory file after run.
      */
      void output() override;

   protected:

      /**
      * Write data that should appear once, at beginning of the file.
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

      using AnalyzerT = typename T::Analyzer;
      using Analyzer<D>::simulator;
      using Analyzer<D>::system;

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;

      /// Has readParam been called?
      long isInitialized_;

   };

}
}
#endif

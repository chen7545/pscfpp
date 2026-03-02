#ifndef RPC_RGRID_TRAJECTORY_READER_H
#define RPC_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"           // base class
#include <prdc/cpu/RField.h>            // member
#include <pscf/math/IntVec.h>           // member
#include <util/containers/DArray.h>     // member
#include <fstream>                      // member
#include <string>

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Trajectory file reader.
   *
   * \ingroup Rpc_Fts_Trajectory_Module
   */
   template <int D>
   class RGridTrajectoryReader : public TrajectoryReader<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      RGridTrajectoryReader<D>(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~RGridTrajectoryReader(){};

      /**
      * Open trajectory file and read header, if any.
      *
      * By convention, this function treats the trajectory filename
      * as the name of an input file, and opens the file using the
      * FileMaster:openInutFile function. This function prepends the
      * input prefix (if any) to the file path.
      *
      * \param filename  trajectory input file name
      */
      void open(std::string filename) override;

      /**
      * Read header of trajectory file (if any).
      */
      void readHeader() override;

      /**
      * Read a single frame from the trajectory file.
      *
      * \return true if a frame is avaiable, false if at end of file
      */
      bool readFrame() override;

      /**
      * Close the trajectory file.
      */
      void close() override;

   protected:

      using TrajectoryReaderT = TrajectoryReader<D>;
      using TrajectoryReaderT::system;

   private:

      // Field configuration
      DArray< RField<D> > wField_;

      // Dimensions of computational mesh (# of points in each direction)
      IntVec<D> meshDimensions_;

      // Trajectory file.
      std::ifstream inputfile_;

      // Has wField_ been allocated?
      bool isAllocated_;

      /**
      * Allocate required memory.
      */
      void allocate();

   };

   // Explicit instantiation declarations
   extern template class RGridTrajectoryReader<1>;
   extern template class RGridTrajectoryReader<2>;
   extern template class RGridTrajectoryReader<3>;

}
}
#endif

#ifndef RPG_STEP_LOGGER_H
#define RPG_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"       // base class template

namespace Pscf {
namespace Rpg {

   template <int D> class Simulator;
   template <int D> class System;

   using namespace Util;

   /**
   * Periodically write snapshots to a trajectory file.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class StepLogger : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      StepLogger(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~StepLogger()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep step index
      */
      virtual void sample(long iStep);

   };

   // Explicit instantiation declarations
   extern template class StepLogger<1>;
   extern template class StepLogger<2>;
   extern template class StepLogger<3>;

}
}
#endif

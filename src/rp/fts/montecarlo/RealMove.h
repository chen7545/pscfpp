#ifndef RP_REAL_MOVE_H
#define RP_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <util/containers/DArray.h>          // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * RealMove generates spatially uncorrelated random field changes.
   *
   * \see \ref rp_RealMove_page "Manual Page"
   *
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class RealMove : public T::McMove
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      RealMove(typename T::McSimulator& simulator);

      /**
      * Destructor.
      */
      ~RealMove();

      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before the simulation loop.
      */
      void setup() override;

      /**
      * Output time contributions.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) override;

   protected:

      using McMoveT = typename T::McMove;
      using RFieldT = typename T::RField;
      using McMoveT::system;
      using McMoveT::simulator;
      using McMoveT::vecRandom;

      /**
      * Attempt unconstrained move.
      */
      void attemptMove();

   private:

      /// New field values, indexed by monomer type.
      DArray< RieldT > w_;

      /// Change in one field component.
      RieldT dwc_;

      /// Standard deviation of field changes.
      double sigma_;

      /// Has memory been allocated?
      bool isAllocated_;

   };

}
}
#endif

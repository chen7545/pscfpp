#ifndef RPG_REAL_MOVE_H
#define RPG_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <prdc/cuda/RField.h>                // member
#include <util/containers/DArray.h>          // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * RealMove is a Monte Carlo move in real space
   *
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class RealMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      RealMove(McSimulator<D>& simulator);

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
      */
      void outputTimers(std::ostream& out) override;

   protected:

      using McMove<D>::system;
      using McMove<D>::simulator;
      using McMove<D>::vecRandom;

      /**
      *  Attempt unconstrained move.
      */
      void attemptMove();

   private:

      /// New field values.
      DArray< RField<D> > w_;

      /// Change in one component of wc.
      RField<D> dwc_;

      // Standard deviation of field changes.
      double sigma_;

      // Has memory been allocated?
      bool isAllocated_;

   };

   // Explicit instantiation declarations
   extern template class RealMove<1>;
   extern template class RealMove<2>;
   extern template class RealMove<3>;

}
}
#endif

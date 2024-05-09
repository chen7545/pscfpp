#ifndef RPC_WC_MOVE_H
#define RPC_WC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <prdc/cpu/RField.h>                 // member
#include <util/containers/DArray.h>          // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * WcMove is a Monte Carlo move in real space 
   *
   * \ingroup Rpc_Simulate_McMove_Module
   */
   template <int D>
   class WcMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      WcMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~WcMove();

      /**
      * Read required parameters from file.
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);
      
      /**
      * Output statistics for this move (at the end of simulation)
      */
      void output();
      
      /**
      * Setup before the beginning of each simulation run
      */
      void setup();
      
      /**
      * Return real move times contributions.
      */
      void outputTimers(std::ostream& out);
      
      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
      using McMove<D>::clearTimers;
      using ParamComposite::read;
      using ParamComposite::setClassName;

   protected:
      
      using McMove<D>::system;
      using McMove<D>::random;
      using McMove<D>::simulator;

      /**
      *  Attempt unconstrained move.
      *
      *  This function should modify the system w fields in r-grid
      *  format, as returned by system().w().rgrid(), in order apply
      *  an unconstrained attempted move. The compressor will then be
      *  applied in order to restore the density constraint.
      *
      */
      void attemptMove();

   private:
      
      /// Move step size, step is selected from [-stepSize_, stepSize_]
      double stepSize_;
      
      /// Change in one component of wc
      RField<D> dwc_;
      
      /// wField after attempted McMove. used in attemptMove() function
      DArray< RField<D> > w_;
      
      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
   
   };
      
   #ifndef RPC_WC_MOVE_TPP
   // Suppress implicit instantiation
   extern template class WcMove<1>;
   extern template class WcMove<2>;
   extern template class WcMove<3>;
   #endif

}
}
#endif

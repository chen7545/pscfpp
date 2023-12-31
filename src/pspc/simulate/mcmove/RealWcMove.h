#ifndef PSPC_REAL_WC_MOVE_H
#define PSPC_REAL_WC_MOVE_H

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
namespace Pspc {

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * RealWcMove is a Monte Carlo move in real space and 
   * change in eigenvector component w fields
   *
   * \ingroup Pspc_Simulate_McMove_Module
   */
   template <int D>
   class RealWcMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      RealWcMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~RealWcMove();

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
      using McMove<D>::mcSimulator;
      using McMove<D>::readProbability;
      using McMove<D>::clearTimers;
      using ParamComposite::read;
      using ParamComposite::setClassName;

   protected:
      
      using McMove<D>::system;
      using McMove<D>::random;

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
      
      // Move step size, step is selected from [-stepSize_, stepSize_]
      double stepSize_;
      
      // Change in one component of wc
      RField<D> dwc_;
      
      // Local copy of w fields
      DArray< RField<D> > w_;
      
      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
   
   };
      
   #ifndef PSPC_REAL_WC_MOVE_TPP
   // Suppress implicit instantiation
   extern template class RealWcMove<1>;
   extern template class RealWcMove<2>;
   extern template class RealWcMove<3>;
   #endif

}
}
#endif

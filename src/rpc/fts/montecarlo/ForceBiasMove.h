#ifndef RPC_FORCE_BIAS_MOVE_H
#define RPC_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
   * ForceBiasMove attempts a Brownian dynamics move.
   *
   * This class implements a Monte Carlo move in which the unconstrained
   * attempted move is created by an explicit Euler Brownian dynamics
   * step. 
   * 
   * Because the probability of attempting a move is not equal to that
   * of generating the reverse move, the acceptance criterion used in
   * the move() function must take into account the ratio of generation 
   * probabilities.
   *
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rpc_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      ForceBiasMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~ForceBiasMove();

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
      * Attempt and accept or reject force bias Monte-Carlo move.
      *
      * \return true if accepted, false if rejected
      */
      bool move();

      /**
      * Return real move times contributions.
      */
      void outputTimers(std::ostream& out);
      
      /**
      * Decide whether dc fields need to be saved for move
      */
      bool needsDc();

   protected:

      // Alias for McMove base class
      using McMoveT = McMove<D>;

      // Protected inherited member functions
      using McMoveT::system;
      using McMoveT::simulator;

   private:

      /// Local copy of w fields
      DArray< RField<D> > w_;

      /// Copy of initial dc field 
      DArray< RField<D> > dc_;

      /// Change in wc
      DArray<RField<D> >  dwc_;

      /// Normal-distributed random field
      RField<D> eta_;

      /// Field used to compute bias for Metropolis criterion
      RField<D> biasField_;

      /// Prefactor of -dc_ in deterministic drift term
      double mobility_;

      /**
      * Compute force bias
      */
      void computeForceBias(RField<D>& result, 
                            RField<D> const & di, 
                            RField<D> const & df, 
                            RField<D> const & dwc, 
                            double mobility);

   };
   
   // Public inline methods

   /*
   * Return whether dc fields need to be saved for ForceBiasMove.
   */
   template <int D>
   inline bool ForceBiasMove<D>::needsDc()
   {  return true; }

   // Explicit instantiation declarations
   extern template class ForceBiasMove<1>;
   extern template class ForceBiasMove<2>;
   extern template class ForceBiasMove<3>;

}
}
#endif

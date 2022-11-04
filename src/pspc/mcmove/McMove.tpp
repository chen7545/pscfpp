#ifndef PSPC_MC_MOVE_TPP
#define PSPC_MC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"

#include <pspc/System.h>
#include <util/archives/Serializable_includes.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMove<D>::McMove(McSimulator<D>& simulator) 
    : simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      randomPtr_(&(systemPtr_->random())),
      nAttempt_(0),
      nAccept_(0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   McMove<D>::~McMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void McMove<D>::readParameters(std::istream &in)
   {}
   
   /*
   * Read the probability from file.
   */
   template <int D>
   void McMove<D>::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }
  
   /*
   * Setup at beginning of loop.
   *
   * Trivial default implementation - initializes counters.
   */
   template <int D>
   void McMove<D>::setup()
   { 
      nAttempt_ = 0;
      nAccept_  = 0;
   }

   /*
   * Trivial default implementation - always returns false.
   */
   template <int D>
   bool McMove<D>::move()
   { 
      incrementNAttempt();

      // Get current Hamiltonian
      double oldHamiltonian = simulator().mcHamiltonian();

      // Save current state 
      simulator().saveMcState();

      // Attempt modification
      attemptMove();

      // Call compressor
      // system().compressor().compress(system.w()));

      // Evaluate new Hamiltonian
      simulator().computeWC();
      simulator().computeMcHamiltonian();
      double newHamiltonian = simulator().mcHamiltonian();

      // Accept or reject move
      bool accept = false;
      double weight = exp(-(newHamiltonian - oldHamiltonian));
      accept = random().metropolis(weight);
      if (accept) {
          incrementNAccept();
      } else {
          simulator().restoreMcState();
      }

      return accept;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void McMove<D>::output()
   {}

}
}
#endif

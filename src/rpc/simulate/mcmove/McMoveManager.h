#ifndef RPC_MC_MOVE_MANAGER_H
#define RPC_MC_MOVE_MANAGER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                      // base class template parameter
#include <util/param/Manager.h>          // base class template
#include <util/containers/DArray.h>      // member template

namespace Util { class Random; }

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class System;
   template <int D> class McSimulator;

   /**
   * Manager for a set of McMove objects.
   *
   * \ingroup Rpc_Simulate_McMove_Module
   */
   template <int D>
   class McMoveManager : public Manager< McMove<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      * \param system parent System
      */
      McMoveManager(McSimulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      ~McMoveManager();

      /**
      * Read instructions for creating McMove objects.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Initialize at beginning of system run.
      *
      * This method calls the initialize method for every McMove.
      */
      void setup();

      /**
      * Choose an McMove at random, using specified probabilities.
      *
      * \return chosen McMove
      */
      McMove<D>& chooseMove();

      /**
      * Return the most recently chosen McMove by const reference.
      *
      * Returns the McMove<D> chosen by the most recent call to
      * chooseMove. Throws an Exception if no move has been chosen.
      *
      * \return previously chosen McMove
      */
      McMove<D> const & chosenMove() const;

      /**
      * Output statistics for all moves.
      */
      void output() const;

      /**
      * Return probability of move i.
      *
      * \param i index for McMove
      * \return probability of McMove number i
      */
      double probability(int i) const;
      
      using Manager< McMove<D> >::size;
      
      /**
      * Log output timing results 
      */
      void outputTimers(std::ostream& out) const;
      
      /**
      * Clear timers 
      */
      void clearTimers();
      
      /**
      * Decide whether any move needs to store cc fields.
      */
      bool needsCc();
      
      /**
      * Decide whether any move needs to store dc fields.
      */
      bool needsDc();
 
   protected:
      
      using Manager< McMove<D> >::setClassName;
   
   private:

      // Private data members

      /**
      * Array of McMove probabilities.
      */
      DArray<double>  probabilities_;
      
      /**
       * Pointer to parent Simulator
       */
      McSimulator<D>* simulatorPtr_;

      /**
      * Pointer to parent System.
      */
      System<D>* systemPtr_;

      /**
      * Pointer to random number generator.
      */
      Random* randomPtr_;

      /**
      * Pointer to most recent McMove.
      */
      McMove<D>* movePtr_;

      // Private member functions

      /**
      * Return pointer to a new McMoveFactory.
      */
      virtual Factory< McMove<D> >* newDefaultFactory() const;

   };

   // Inline functions

   /*
   * Return probability of move number i
   */
   template <int D>
   inline double McMoveManager<D>::probability(int i) const
   {
      assert(i >= 0);  
      assert(i < size());  
      return probabilities_[i];
   }

   #ifndef RPC_MC_MOVE_MANAGER_TPP
   // Suppress implicit instantiation
   extern template class McMoveManager<1>;
   extern template class McMoveManager<2>;
   extern template class McMoveManager<3>;
   #endif

}
}
#endif

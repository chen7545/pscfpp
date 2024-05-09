#ifndef RPC_WC_MOVE_TPP
#define RPC_WC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WcMove.h"
#include "McMove.h" 
#include <rpc/simulate/mcmove/McSimulator.h>
#include <rpc/System.h>
#include <util/param/ParamComposite.h>
#include <util/archives/Serializable_includes.h>
#include <util/random/Random.h>


namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   WcMove<D>::WcMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      isAllocated_(false)
   { setClassName("WcMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   WcMove<D>::~WcMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void WcMove<D>::readParameters(std::istream &in)
   {
      // Read the probability
      readProbability(in);
      
      // Attampt move range [A, -A]
      read(in, "A", stepSize_);
      
      // Allocate memory for private containers
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      if (!isAllocated_){
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(meshDimensions);
         }
         dwc_.allocate(meshDimensions);
         isAllocated_ = true;
      }
      
   }
   
   template <int D>
   void WcMove<D>::setup()
   {  
      McMove<D>::setup();
      
      // Check array capacities
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
   }
   
   
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void WcMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // Copy current W fields from parent system into w_
      for (int i = 0; i < nMonomer; ++i) {
         w_[i] = system().w().rgrid(i);
      }
      
      double evec;
      // eigenvectors of projected chi matrix
      for (int j = 0; j < nMonomer - 1; j++){
         for (int k = 0; k < meshSize; k++){
            dwc_[k] = random().uniform(-stepSize_,stepSize_);
         }
         
         // Loop over monomer types (back to monomer fields)
         for (int i = 0; i < nMonomer; i++){
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            for (int k = 0; k < meshSize; ++k) {
               w[k] += evec*dwc_[k];
            }
         }
      }
      
      // Set modified fields in parent system
      system().setWRGrid(w_);
   }


   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void WcMove<D>::output()
   {}
   
   template<int D>
   void WcMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Wc Move times contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif

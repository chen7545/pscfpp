#ifndef PSPC_REAL_WC_MOVE_TPP
#define PSPC_REAL_WC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealWcMove.h"
#include "McMove.h" 
#include <pspc/simulate/mcmove/McSimulator.h>
#include <util/param/ParamComposite.h>
#include <pspc/System.h>
#include <util/archives/Serializable_includes.h>
#include <util/random/Random.h>


namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealWcMove<D>::RealWcMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      isAllocated_(false)
   { setClassName("RealWcMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   RealWcMove<D>::~RealWcMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void RealWcMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", stepSize_);
   }
   
   template <int D>
   void RealWcMove<D>::setup()
   {  
      McMove<D>::setup();
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
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void RealWcMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
       // Copy current W fields from parent system into wc_
      for (int i = 0; i < nMonomer; i++) {
         w_[i] = system().w().rgrid(i);
      }
      
      // Evaluate component property
      mcSimulator().computeWc();
      
      // For multi-component copolymer
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (int j = 0; j < nMonomer - 1; j++){
         for (int k = 0; k < meshSize; k++){
            //Random number generator
            double r = random().uniform(-stepSize_,stepSize_);
            dwc_[k] = r;
         }
         // Loop over monomer types to modify monomer w field
         for (int i = 0; i < nMonomer; i++) {
            RField<D> & wn = w_[i];
            evec = mcSimulator().chiEvecs(j,i);
            for (int k = 0; k < meshSize; k++){
               wn[k] += evec * dwc_[k];
            }
         }
      }
      // Set modified fields in parent system
      system().setWRGrid(w_);
      mcSimulator().clearData();
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void RealWcMove<D>::output()
   {}
   
   template<int D>
   void RealWcMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "RealWc Move times contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif

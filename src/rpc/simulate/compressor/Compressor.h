#ifndef RPC_COMPRESSOR_H
#define RPC_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>
#include <iostream>
#include <vector>

namespace Pscf {
namespace Rpc
{

   template <int D> class System;

   using namespace Util;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Rpc_Simulate_Compressor_Module
   */
   template <int D>
   class Compressor : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      Compressor(System<D>& system);

      /**
      * Destructor.
      */
      ~Compressor();

      /**
      * Iterate Langrange multiplier field.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int compress() = 0;

      /**
      * Output report of timing results to stream.
      *
      * \param out output stream for results
      */
      virtual void outputTimers(std::ostream& out) = 0;

      /**
      * Clear timers.
      */
      virtual void clearTimers() = 0;
      
      /**
      * Get the ratio of error reduction by AM step 1
      */
      virtual std::vector<double> stepOneRatioVector() = 0;
      
      /**
      * Get the predicted of error reduction by AM step 1
      */
      virtual std::vector<double> predictRatioVector() = 0;
      
      /**
      * Get the ratio of error reduction by AM step 1
      */
      virtual std::vector<double> stepTwoRatioVector() = 0;
      
      /**
      * Return error after sampling step, error at compressor Itr0.
      */
      virtual double errorItr0() = 0;

      /**
      * Get the number of times the MDE has been solved.
      */
      int mdeCounter();
      
      /**
      * Get the total number of iterations 
      */
      int totalItr();

   protected:

      /**
      * Return const reference to parent system.
      */
      System<D> const & system() const
      {  
         UTIL_ASSERT(sysPtr_);
         return *sysPtr_;
      }

      /**
      * Return non-const reference to parent system.
      */
      System<D>& system()
      {  
         UTIL_ASSERT(sysPtr_);
         return *sysPtr_;
      }

      /**
      * Count how many times MDE has been solved.
      */
      int mdeCounter_;
      
      /**
      * Count the total of iteration required to converge
      */
      int totalItr_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Member functions

   // Constructor
   template <int D>
   Compressor<D>::Compressor(System<D>& system)
    : mdeCounter_(0),
      sysPtr_(&system)
   {  setClassName("Compressor"); }

   // Destructor
   template <int D>
   Compressor<D>::~Compressor()
   {}

   // Get number of times MDE has been solved.
   template <int D>
   inline int Compressor<D>::mdeCounter()
   {  return mdeCounter_; }
   
   // Get number of times MDE has been solved.
   template <int D>
   inline int Compressor<D>::totalItr()
   {  return totalItr_; }


} // namespace Rpc
} // namespace Pscf
#endif

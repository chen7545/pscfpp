#ifndef PSPG_COMPRESSOR_H
#define PSPG_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Pspg_Compressor_Module
   */
   template <int D>
   class Compressor : public ParamComposite
   {

   public:
      #if 0
      /**
      * Default constructor.
      */
      Compressor();
      #endif

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
      
      virtual int counterMDE() = 0;
      
      /**
      * Count how many times error goes up
      */
      virtual int counterErrorUp() = 0;
      
      /**
      * Log output timing results 
      */
      virtual void outputTimers(std::ostream& out) = 0;
      
      /**
      * Clear timers 
      */
      virtual void clearTimers() = 0;

      /**
      * Return const reference to parent system.
      */
      System<D> const & system() const
      {  return *sysPtr_;}

   protected:

      /**
      * Return reference to parent system.
      */
      System<D>& system() 
      {  return *sysPtr_;}

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Inline member functions
   #if 0
   // Default constructor
   template <int D>
   inline Compressor<D>::Compressor()
    : sysPtr_(&system)
   {  setClassName("Compressor"); }
   #endif

   // Constructor
   template <int D>
   Compressor<D>::Compressor(System<D>& system)
    : sysPtr_(&system)
   {  setClassName("Compressor"); }

   // Destructor
   template <int D>
   Compressor<D>::~Compressor()
   {}
   
} // namespace Pspg
} // namespace Pscf
#endif

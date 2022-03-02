#ifndef PSPG_ITERATOR_FACTORY_TPP
#define PSPG_ITERATOR_FACTORY_TPP

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include "AmIterator.h"
#include "AmIteratorOld.h"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D>::IteratorFactory(System<D>& system)
    : sys_(&system)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<D>* IteratorFactory<D>::factory(const std::string &className) const
   {
      Iterator<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "AmIterator") {
         ptr = new AmIterator<D>(*sys_);
      } else if (className == "AmIteratorOld") {
         ptr = new AmIteratorOld<D>(*sys_);
      }

      return ptr;
   }

}
}
#endif
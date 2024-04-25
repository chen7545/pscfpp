#ifndef RPC_COMPRESSOR_FACTORY_TPP
#define RPC_COMPRESSOR_FACTORY_TPP

#include "CompressorFactory.h"  

// Subclasses of Compressor 
#include "AmCompressor.h"
#include "LrCompressor.h"
#include "LrAmCompressor.h"
#include "LrPostAmCompressor.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   CompressorFactory<D>::CompressorFactory(System<D>& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Compressor subclass className.
   */
   template <int D>
   Compressor<D>* 
   CompressorFactory<D>::factory(std::string const &className) const
   {
      Compressor<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "Compressor" || className == "AmCompressor") {
         ptr = new AmCompressor<D>(*sysPtr_);
      } else if (className == "LrCompressor") {
         ptr = new LrCompressor<D>(*sysPtr_);
      } else if (className == "LrAmCompressor") {
         ptr = new LrAmCompressor<D>(*sysPtr_);
      } else if (className == "LrPostAmCompressor") {
         ptr = new LrPostAmCompressor<D>(*sysPtr_);
      }
      
      return ptr;
   }

}
}
#endif

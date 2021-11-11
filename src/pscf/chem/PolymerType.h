#ifndef PSCF_POLYMER_TYPE_H
#define PSCF_POLYMER_TYPE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#include <iostream>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Struct containing an enumeration of polymer structure types.
   */
   struct PolymerType {
      enum Enum {Branched, Linear};
   };

   /**
   * Input stream extractor for a PolymerType enumeration.
   *
   * \param in input stream
   * \param type value of PolymerType to be read from file
   */ 
   std::istream& operator >> (std::istream& in, PolymerType::Enum& type); 

   /**
   * Input stream extractor for a PolymerType enumeration.
   *
   * \param out  output stream
   * \param type  value of PolymerType to be written 
   */ 
   std::ostream& operator << (std::ostream& out, PolymerType::Enum& type); 

   /**
   * Input stream extractor for a PolymerType enumeration.
   *
   * \param out  output stream
   * \param type  value of PolymerType to be written 
   */ 
   template <class Archive>
   inline void 
   serialize(Archive& ar, PolymerType::Enum& data, const unsigned int version)
   {  serializeEnum(ar, data, version); }

}
#endif
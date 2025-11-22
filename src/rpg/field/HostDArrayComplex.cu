/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostDArrayComplex.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Default constructor.
   */
   HostDArrayComplex::HostDArrayComplex() 
    : HostDArray<ValueType>()
   {}

   /*
   * Allocating constructor.
   */
   HostDArrayComplex::HostDArrayComplex(int capacity)
    : HostDArray<ValueType>(capacity)
   {}

   /*
   * Copy constructor.
   */
   HostDArrayComplex::HostDArrayComplex(HostDArrayComplex const& other)
    : HostDArray<ValueType>(other)
   {}

   /*
   * Destructor.
   */
   HostDArrayComplex::~HostDArrayComplex()
   {}

} // namespace Rpg
} // namespace Pscf

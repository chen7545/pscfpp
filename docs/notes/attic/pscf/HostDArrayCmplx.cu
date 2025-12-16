/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostDArrayCmplx.h"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   /*
   * Default constructor.
   */ 
   HostDArrayCmplx()
    : Base()
   {}

   /*
   * Allocating constructor.
   */ 
   HostDArrayCmplx(int n)
    : Base(n)
   {}

   /*
   * Destructor.
   */ 
   ~HostDArrayCmplx()
   {}

} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif

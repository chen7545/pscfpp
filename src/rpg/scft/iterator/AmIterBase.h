#ifndef RPG_AM_ITER_BASE_H
#define RPG_AM_ITER_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>  // base class template
#include <rpg/scft/iterator/Iterator.h>    // base class argument
#include <pscf/cuda/DeviceArray.h>         // base class argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<1>, DeviceArray<cudaReal> >;
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<2>, DeviceArray<cudaReal> >;
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<3>, DeviceArray<cudaReal> >;
}
#endif

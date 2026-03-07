#ifndef RPG_AM_ITERATOR_DEV_H
#define RPG_AM_ITERATOR_DEV_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* Declarations of explicit instantiations of AmiteratorTmpl that 
* used as base classes for class template Rpg::AmIteratorGrid.
*/

#include <pscf/iterator/AmIteratorTmpl.h>  // base class template
#include <rpg/scft/iterator/Iterator.h>    // base class argument
#include <pscf/cuda/DeviceArray.h>         // base class argument
#include <pscf/cuda/cudaTypes.h>           // base class argument

namespace Pscf {
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<1>, DeviceArray<cudaReal> >;
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<2>, DeviceArray<cudaReal> >;
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<3>, DeviceArray<cudaReal> >;

}
#endif

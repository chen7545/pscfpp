#ifndef RPG_AM_ITERATOR_HOST_H
#define RPG_AM_ITERATOR_HOST_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* Declarations of explicit instantiations of AmIteratorTmpl that are
* used as base classes for class template Rpg::AmIteratorBasis.
*/

#include <pscf/iterator/AmIteratorTmpl.h>  // base class template
#include <rpg/scft/iterator/Iterator.h>    // base class argument
#include <util/containers/DArray.h>        // base class argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<1>, DArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<2>, DArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpg::Iterator<3>, DArray<double> >;
}
#endif

#ifndef RPC_AM_ITER_BASE_H
#define RPC_AM_ITER_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>   // base class template
#include <rpc/scft/iterator/Iterator.h>     // base class argument
#include <util/containers/DArray.h>         // base class argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class AmIteratorTmpl< Rpc::Iterator<1>, DArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpc::Iterator<2>, DArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpc::Iterator<3>, DArray<double> >;
}
#endif

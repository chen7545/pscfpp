/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorHost.h"
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {

   template 
   class AmIteratorTmpl< Rpg::Iterator<1>, DArray<double> >;
   template 
   class AmIteratorTmpl< Rpg::Iterator<2>, DArray<double> >; 
   template 
   class AmIteratorTmpl< Rpg::Iterator<3>, DArray<double> >;

}

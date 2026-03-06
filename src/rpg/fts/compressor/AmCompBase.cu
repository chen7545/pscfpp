/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompBase.h"
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template 
   class AmIteratorTmpl< Rpc::Compressor<1>, DeviceArray<cudaReal> >;
   template 
   class AmIteratorTmpl< Rpc::Compressor<2>, DeviceArray<cudaReal> >; 
   template 
   class AmIteratorTmpl< Rpc::Compressor<3>, DeviceArray<cudaReal> >;
}

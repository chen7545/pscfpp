#ifndef RPG_AM_COMPRESSOR_H
#define RPG_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/AmCompressor.h> // direct base class template
#include <rpg/system/Types.h>               // direct base template argument
#include <prdc/cuda/RField.h>               // direct base class member
#include <pscf/iterator/AmIteratorTmpl.h>   // indirect base class template
#include <rpg/fts/compressor/Compressor.h>  // indirect base argument
#include <pscf/cuda/DeviceArray.h>          // indirect base argument
#include <pscf/cuda/cudaTypes.h>            // indirect base argument

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Anderson mixing compressor.
   *
   * \see \ref rp_AmCompressor_page "Manual Page"
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor
    : public Rp::AmCompressor<D, Rpg::Types<D>, DeviceArray<cudaReal> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      AmCompressor(System<D>& system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class AmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      extern template 
      class AmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      extern template 
      class AmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
   namespace Rpg {
      extern template class AmCompressor<1>;
      extern template class AmCompressor<2>;
      extern template class AmCompressor<3>;
   }
}
#endif

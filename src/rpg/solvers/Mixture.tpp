#ifndef RPG_MIXTURE_TPP
#define RPG_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"               // class header

#include "Polymer.h"
#include "Solvent.h"
#include "Block.h"
#include "Propagator.h"
#include <rpg/field/FieldIo.h>
#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/VecOp.h>

#include <rp/solvers/Mixture.tpp>  // base class template implementation

namespace Pscf {
namespace Rpg {

   using namespace Prdc;

   /*
   * Constructor
   */
   template <int D>
   Mixture<D>::Mixture()
    : Rp::Mixture<D, Types<D> >(),
      useBatchedFFT_(true)
   {}

   /*
   * Read all parameters and initialize.
   */
   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      RpMixtureT::readParameters(in);

      // Optionally read useBatchedFFT boolean
      useBatchedFFT_ = true;
      ParamComposite::readOptional(in, "useBatchedFFT", useBatchedFFT_);
   }

   /*
   * Allocate memory for all blocks.
   */
   template <int D>
   void Mixture<D>::allocateBlocks()
   {
      const double ds = RpMixtureT::ds();
      const int np = MixtureBase<cudaReal>::nPolymer();
      int i, j;
      for (i = 0; i < np; ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).allocate(ds, useBatchedFFT_);
         }
      }
   }

} // namespace Rpg
} // namespace Pscf
#endif

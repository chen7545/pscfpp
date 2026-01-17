#ifndef RPG_SIMULATOR_TPP
#define RPG_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>
#include <rpg/fts/simulator/SimState.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/fts/compressor/CompressorFactory.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/fts/perturbation/PerturbationFactory.h>
#include <rpg/fts/ramp/Ramp.h>
#include <rpg/fts/ramp/RampFactory.h>

#include <prdc/cuda/RField.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/math/IntVec.h>

#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>

#include <gsl/gsl_eigen.h>

#include <rp/fts/simulator/Simulator.tpp>  // base class implementation

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Simulator<D>::Simulator(typename Types<D>::System& system)
    : Base(system, *this)
   {}

   template <int D>
   void Simulator<D>::initializeVecRandom()
   {  Base::vecRandom().setSeed(Base::seed_); }

}
}
#endif

#ifndef RPC_SIMULATOR_TPP
#define RPC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>
#include <rpc/fts/simulator/SimState.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/fts/compressor/CompressorFactory.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/perturbation/PerturbationFactory.h>
#include <rpc/fts/ramp/Ramp.h>
#include <rpc/fts/ramp/RampFactory.h>

#include <prdc/cpu/RField.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/math/IntVec.h>

#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>
#include <gsl/gsl_eigen.h>

#include <rp/fts/simulator/Simulator.tpp>  // base class implementation

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Simulator<D>::Simulator(typename Types<D>::System& system)
    : Base(system, *this)
   {  Base::vecRandom().associate(Base::random()); }

}
}
#endif

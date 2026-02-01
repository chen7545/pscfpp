#ifndef RPC_MC_SIMULATOR_TPP
#define RPC_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <rpc/fts/montecarlo/McMove.h>
#include <rpc/fts/montecarlo/McMoveFactory.h>
#include <rpc/fts/montecarlo/McMoveManager.h>
#include <rpc/fts/analyzer/Analyzer.h>
#include <rpc/fts/analyzer/AnalyzerFactory.h>
#include <rpc/fts/analyzer/AnalyzerManager.h>
#include <rpc/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpc/fts/trajectory/TrajectoryReader.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/ramp/Ramp.h>
#include <rpc/system/System.h>
#include <prdc/cpu/RField.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/montecarlo/McSimulator.tpp>

namespace Pscf {
namespace Rpc {

   /*
   * Constructor.
   */
   template <int D>
   McSimulator<D>::McSimulator(System<D>& system)
    : Rp::McSimulator<D, Types<D> >(system, *this)
   {}

}
}
#endif

#ifndef RPG_MC_SIMULATOR_TPP
#define RPG_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <rpg/fts/montecarlo/McMove.h>
#include <rpg/fts/montecarlo/McMoveFactory.h>
#include <rpg/fts/montecarlo/McMoveManager.h>
#include <rpg/fts/analyzer/Analyzer.h>
#include <rpg/fts/analyzer/AnalyzerFactory.h>
#include <rpg/fts/analyzer/AnalyzerManager.h>
#include <rpg/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpg/fts/trajectory/TrajectoryReader.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/fts/ramp/Ramp.h>
#include <rpg/system/System.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/CudaVecRandom.h>

#include <rp/fts/montecarlo/McSimulator.tpp>

namespace Pscf {
namespace Rpg {

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

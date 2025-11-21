#ifndef RPC_TYPES_H
#define RPC_TYPES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declarations
namespace Pscf {
   class Interaction;
   template <typename WT> class Species;
   template <typename WT> class PolymerSpecies;
   template <typename WT> class SolventSpecies;
   template <typename WT> class MixtureBase;
   namespace Prdc {
      class Environment;
      namespace Cpu {
         template <int D> class CField;
         template <int D> class FFT;
         template <int D> class CFieldComparison;
         template <int D> class WaveList;
      }
   }
   namespace Rpc {
      template <int D> class System;
      template <int D> class Mixture;
      // template <int D> class MixtureModifier;
      template <int D> class Polymer;
      template <int D> class Solvent;
      template <int D> class Block;
      template <int D> class Propagator;
      template <int D> class Domain;
      template <int D> class FieldIo;
      template <int D> class WFields;
      template <int D> class CFields;
      // template <int D> class Simulator;
      // template <int D> class SimulatorFactory;
   }
}

namespace Pscf {
namespace Rpc {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * List of aliases for types used the in Rpc namespace.
   *
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class Types
   {

   public:

      using System = Rpc::System<D>;

      using Mixture = Rpc::Mixture<D>;
      // using MixtureModifier = Rpc::MixtureModifier<D>;
      using Polymer = Rpc::Polymer<D>;
      using Solvent = Rpc::Solvent<D>;
      using Block = Rpc::Block<D>;
      using Propagator = Rpc::Propagator<D>;

      using Species = Pscf::Species<double>;
      using PolymerSpecies = Pscf::PolymerSpecies<double>;
      using SolventSpecies = Pscf::SolventSpecies<double>;
      using MixtureBase = Pscf::MixtureBase<double>;

      using Interaction = Pscf::Interaction;
      using Domain = Rpc::Domain<D>;
      using FieldIo = Rpc::FieldIo<D>;
      using ScftThermo = Rpc::ScftThermo<D>;

      using WFields = Rpc::WFields<D>;
      using CFields = Rpc::CFields<D>;
      using Mask = Rpc::Mask<D>;

      // using Simulator = Rpc::Simulator<D>;
      // using SimulatorFactory = Rpc::SimulatorFactory<D>;

      using CField = Prdc::Cpu::CField<D>;
      using FFT = Prdc::Cpu::FFT<D>;
      using CFieldComparison = Prdc::Cpu::CFieldComparison<D>;
      using WaveList = Prdc::Cpu::WaveList<D>;

   };

} // namespace Rpc
} // namespace Pscf
#endif

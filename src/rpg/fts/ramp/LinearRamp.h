#ifndef RPG_LINEAR_RAMP_H
#define RPG_LINEAR_RAMP_H

#include <rpg/fts/ramp/Ramp.h>           // base class
#include <rpg/fts/ramp/RampParameter.h>  // member (templ parameter)
#include <util/containers/DArray.h>      // member (template)

namespace Pscf {
namespace Rpg {

   using namespace Util;

   template <int D> class Simulator;

   /**
   * Linear ramp - parameters vary linearly with step index.
   *
   * \see \ref rp_LinearRamp_page "Manual Page"
   * \see \ref psfts_ramp_page "Manual Page"
   * \ingroup Rpg_Fts_Ramp_Module
   */
   template <int D>
   class LinearRamp : public Ramp<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      LinearRamp(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~LinearRamp();

      /**
      * Read parameters from parameter file input stream.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Set nStep and complete initialization.
      *
      * This method is called just before the beginning of the main
      * simulation loop.
      *
      * \param nStep number of steps planned for this simulation
      */
      void setup(int nStep) override;

      /**
      * Set new parameters values in associated System and Simulator.
      *
      * \param iStep  current simulation step index
      */
      void setParameters(int iStep) override;

      /**
      * Output information at the end of a simulation.
      */
      void output() override;

   protected:

      using RampT = Ramp<D>;
      using RampParameterT = RampParameter<D>;
      using Ramp<D>::nStep_;
      using Ramp<D>::simulator;

   private:

      // Number of variable parameters
      int nParameter_;

      // Array of variable parameters
      DArray< RampParameterT > parameters_;

   };

   // Explicit instantiation declarations
   extern template class LinearRamp<1>;
   extern template class LinearRamp<2>;
   extern template class LinearRamp<3>;

}
}
#endif

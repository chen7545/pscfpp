#ifndef RPG_AVERAGE_ANALYZER_H
#define RPG_AVERAGE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                         // grandparent base class
#include <rp/fts/analyzer/AverageAnalyzer.h>  // base class template
#include <rpg/system/Types.h>                 // class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of a single sampled real variables,
   * and optionally writes values or block averages to a data file during a
   * simulation.  It is intended for use as a base class for any Analyzer
   * that computes and evaluates an average for a single physical variable.
   *
   * This class is basically a named instantiation of the base class
   * template Rp::AverageAnalyzer, using aliases defined in Rpg::Types<D>
   * to specialize to types used in the Rpg namespace. See documentation
   * of the base class for details.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class AverageAnalyzer : public Rp::AverageAnalyzer<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      AverageAnalyzer(Simulator<D>& simulator, System<D>& system);

      /// Alias for base class
      using AnalyzerT = typename Types<D>::Analyzer;
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageAnalyzer<1, Rpg::Types<1> >;
      extern template class AverageAnalyzer<2, Rpg::Types<2> >;
      extern template class AverageAnalyzer<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class AverageAnalyzer<1>;
      extern template class AverageAnalyzer<2>;
      extern template class AverageAnalyzer<3>;
   }
}
#endif

#ifndef RPC_WC_FOURIER_MOVE_H
#define RPC_WC_FOURIER_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <prdc/cpu/RField.h>                 // member
#include <prdc/cpu/RFieldDft.h>              // member
#include <util/containers/DArray.h>          // member
//#include <prdc/crystal/shiftToMinimum.h>
//#include <pscf/math/IntVec.h>
//#include <util/param/ParamComposite.h>
//#include <util/global.h>

namespace Pscf {
namespace Rpc
{

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * WcFourierMove is a Monte Carlo move in Fourier space
   *
   * \ingroup Rpc_Simulate_McMove_Module
   */
   template <int D>
   class WcFourierMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      WcFourierMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~WcFourierMove();

      /**
      * Read required parameters from file.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream &in);
      
      /**
      * Output statistics for this move (at the end of simulation)
      */
      void output();
      
      /**
      * Setup before the beginning of each simulation run
      */
      void setup();
      
      /**
      * Return fourier move times contributions.
      *
      * \param out  output stream for timer statistics
      */
      void outputTimers(std::ostream& out);
      
      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
      using McMove<D>::clearTimers;
      using ParamComposite::read;
      using ParamComposite::setClassName;

   protected:
      
      using McMove<D>::system;
      using McMove<D>::random;
      using McMove<D>::simulator;
      
      /**
      * Attempt unconstrained move.
      *
      * This function should convert system w fields from r-grid format to
      * RFieldDft format, modify the system w fields in fourier space, and
      * convert back to r-grid format, as returned by system().w().rgrid(),
      * in order apply an unconstrained attempted move. The compressor will 
      * then be applied in order to restore the density constraint.
      */
      void attemptMove();


   private:
      
      /**
      * Compute square of radius of gyration Rg = (Nb^2/6)
      */
      void computeRgSquare();
      
      /**
      * Compute Leibler's F(x,f) function. 
      * 
      * F(x,f) = g(1,x)/{g(f,x)g(1-f,x) - [g(1,x) - g(f,x) - g(1-f,x)]^2/4}
      *
      * \param x nondimensionalized value of q^2 Rg^2
      */
      double computeF(double x);
      
      /**
      * Compute Debye function g(f,x).
      * 
      * Debye function: g(f,x) = 2[fx + exp(-fx)-1]/x^2
      *
      * \param f fraction of A block
      * \param x nondimensionalized value of q^2 Rg^2
      */
      double computeDebye(double f, double x);
      
      /**
      * Compute Fredrickson-Helfand structure for specific qSquare.
      *
      * For diblock copolymers,Fredrickson-Helfand structure 
      * S(q)/N = 1/(F(x,f) - F* + epsilon ) x = q^2Rg^2. 
      * q is wave vector, Rg is radius of gyration
      */
      double computeS(double qSquare);
      
      /**
      * Compute Fredrickson-Helfand structure factor
      */
      void computeStructureFactor();
      
      /**
      * Input variable, Move step size stepSize_.
      *
      * In Fourier space Delta W(q) is randomly selected from a uniform
      * distribution [-stepSize_*S(q)^(1/2), stepSize_*S(q)^(1/2)].
      */ 
      double stepSize_;
     
      /** 
      * Input variable, F*. 
      *
      * Ex: at f = 0.5, F* = 20.990; at f = 0.45, F* = 21.396.
      */
      double fStar_;
      
      /**
      * Input variable tau_. 
      *
      * User can calculate using Equation (5.2) from reference:
      * "Fluctuation effects in the theory of microphase separation
      * in block copolymers." Fredrickson, G. H., & Helfand, E. 
      * J. Chem. Phys. 87(1), 697-705 (1987).
      */
      double tau_;
      
      /// Input variable, diblock volume fraction.
      double f_;
      
      /// Input variable, degree of polymerization N. 
      double n_;

      /// Input variable, statistical segment length b.
      double b_;
      
      /// Radius of gyration 
      double rgSquare_;
      
      /// Fredrickson-Helfand structure
      DArray<double> structureFactors_; 
      
      /// M - 1 Eigenvector composition fields in in fourier space.
      DArray< RFieldDft<D> > wck_;
      
      /**
      * M-1 Eigenvector composition fields in in real space.
      * after attempt McMove.
      */ 
      DArray< RField<D> > wc_;
      
      /// wField in real space after attempt McMove.
      DArray< RField<D> > w_;
      
      /// Has the variable been allocated?
      bool isAllocated_;
   
   };
      
   #ifndef RPC_WC_FOURIER_MOVE_TPP
   // Suppress implicit instantiation
   extern template class WcFourierMove<1>;
   extern template class WcFourierMove<2>;
   extern template class WcFourierMove<3>;
   #endif

}
}
#endif

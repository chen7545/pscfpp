#ifndef RPC_WC_FOURIER_MOVE_TPP
#define RPC_WC_FOURIER_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WcFourierMove.h"
#include "McMove.h" 
#include <rpc/simulate/mcmove/McSimulator.h>
#include <rpc/System.h>      
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
//#include <pscf/math/RealVec.h>
#include <util/random/Random.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <iostream>
#include <complex>
#include <cmath>

//#include <util/archives/Serializable_includes.h>
//#include <util/format/Int.h>

namespace Pscf {
namespace Rpc 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   WcFourierMove<D>::WcFourierMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      isAllocated_(false)
   { setClassName("WcFourierMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   WcFourierMove<D>::~WcFourierMove()
   {}

   /*
   * ReadParameters
   */
   template <int D>
   void WcFourierMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", stepSize_);
      read(in, "F*", fStar_);
      read(in, "tau", tau_);
      read(in, "N", n_);
      read(in, "volumeFraction", f_);
      read(in, "statisticalSegmentLength", b_);
      
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      if (!isAllocated_){
         w_.allocate(nMonomer);
         wc_.allocate(nMonomer - 1);
         wck_.allocate(nMonomer - 1);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(dimensions);
         }
         
         for (int i = 0; i < nMonomer - 1; ++i){
            wc_[i].allocate(dimensions);
            wck_[i].allocate(dimensions);
         }
         
         isAllocated_ = true;
      }
      structureFactors_.allocate(system().mesh().size());
   }
   
   template <int D>
   void WcFourierMove<D>::setup()
   {  
      McMove<D>::setup();
      computeRgSquare();
      computeStructureFactor();
   }
   
   /**
   * Compute radius of gyration Rg^2 = Nb^2/6
   */
   template <int D>
   void WcFourierMove<D>::computeRgSquare()
   {
      rgSquare_ = n_ * b_ * b_ / 6.0;
   }

   /*
   * Debye function g(f,x)
   */
   template <int D>
   double WcFourierMove<D>::computeDebye(double f, double x)
   {
      double fx = f * x;
      return 2.0 * (fx + exp(-fx) - 1.0) / (x * x);
   }
   
   /*
   * F(x,f) = g(1,x)/{g(f,x)g(1-f,x) - [g(1,x) - g(f,x) - g(1-f,x)]^2/4}
   */
   template <int D>
   double WcFourierMove<D>::computeF(double x)
   {
      //Only for diblock
      if (f_ > 0.5){
         f_ = 1.0 - f_;
      }
      
      double g1 = computeDebye(1.0, x);
      double gf = computeDebye(f_, x);
      double g1f = computeDebye(1.0 - f_, x);
      double denominator = gf * g1f - pow(g1 - gf- g1f, 2.0)/4.0;
      return g1 / denominator;
   }
   
   /*
   * S(q)/N Fredrickson-Helfand structure factor. S(q)/N = 1/(F(x) - F* + tau_)
   */  
   template <int D>
   double WcFourierMove<D>::computeS(double qSquare)
   {
      //Input is square magnitude of reciprocal basis vector
      double x = qSquare * rgSquare_;
      double denominator =  computeF(x) - fStar_ + tau_;
      // Nbar = Rho^(-2)*a^6* N && Rho = 1/vMonomer && vMonomer = 1/sqrt(Nbar)
      const double vMonomer = system().mixture().vMonomer();
      return 1.0 / (vMonomer * denominator);
   }
   
   template <int D>
   void WcFourierMove<D>::computeStructureFactor()
   {
      MeshIterator<D> itr;
      itr.setDimensions(system().mesh().dimensions());
      IntVec<D> G; IntVec<D> Gmin; 
      double qSq; double S_; 
      for (itr.begin(); !itr.atEnd(); ++itr){
         if (itr.rank() == 0){
            structureFactors_[itr.rank()] = 0;
         }
         // Obtain square magnitude of reciprocal basis vector
         if (itr.rank() >= 1){
            G = itr.position();
            Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
            qSq = system().unitCell().ksq(Gmin);
            // Compute Fredrickson-Helfand structure factor
            S_ = computeS(qSq);
            structureFactors_[itr.rank()] = S_;
         }
      }
   }
   
   
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void WcFourierMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // Convert real grid wc to KGrid format
      for (int i = 0; i < nMonomer - 1; ++i) {
         system().fft().forwardTransform(simulator().wc(i), wck_[i]);
      }
      
      /// Move step size in Fourier space Delta_W(q) is randomly selected from 
      /// uniform distribution [-stepSize_*S(q)^(1/2), stepSize_*S(q)^(1/2)]      
      double S_; 
      for (int i = 0; i < nMonomer - 1; i++){
         MeshIterator<D> itr;
         itr.setDimensions(wck_[0].dftDimensions());
         for (itr.begin(); !itr.atEnd(); ++itr){
            // Obtain square magnitude of reciprocal basis vector
            if (itr.rank() >= 1){
               // Compute Fredrickson-Helfand structure factor
               S_ = structureFactors_[itr.rank()];
               // Generate random complex number
               double real = random().uniform(-stepSize_,stepSize_);
               double imag = random().uniform(-stepSize_,stepSize_);
               // Attempt Move in ourier (k-grid) format
               wck_[i][itr.rank()][0] += real * sqrt(S_); 
               wck_[i][itr.rank()][1] += imag * sqrt(S_);
            }
         }
      }

      // Convert Fourier (k-grid) format back to Real grid format
      for (int i = 0; i < nMonomer - 1; ++i) {
         system().fft().inverseTransform(wck_[i], wc_[i]);
      }
      // Convert eigenvector compositions to monomer fields
      double evec;
   
      // Convert pressure field 
      for (int i = 0; i < nMonomer; i++){
         evec = simulator().chiEvecs(nMonomer,i);
         for (int k = 0; k < meshSize; ++k) {
            w_[i][k] = evec * simulator().wc(nMonomer-1)[k];
         }
      }
   
      // Convert M - 1 fluctuating fields
      for (int j = 0; j < nMonomer - 1; j++){
         for (int i = 0; i < nMonomer; i++){
            evec = simulator().chiEvecs(j,i);
            for (int k = 0; k < meshSize; ++k) {
               w_[i][k] += evec*wc_[j][k];
            }
         }
      }
      
      // Update attemptMove
      system().setWRGrid(w_);

   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void WcFourierMove<D>::output()
   {}
   
   template<int D>
   void WcFourierMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "WcFourier Move times contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif

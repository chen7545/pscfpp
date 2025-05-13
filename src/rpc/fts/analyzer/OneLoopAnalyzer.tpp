#ifndef RPC_ONE_LOOP_ANALYZER_TPP
#define RPC_ONE_LOOP_ANALYZER_TPP

#include "OneLoopAnalyzer.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/System.h>

#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/Debye.h>
#include <pscf/chem/EdgeIterator.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <fftw3.h>

#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>

//#include <unordered_map>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   OneLoopAnalyzer<D>::OneLoopAnalyzer(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      kSize_(1),
      totalN_(0),
      isInitialized_(false)
   {  setClassName("OneLoopAnalyzer"); }


   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void OneLoopAnalyzer<D>::readParameters(std::istream& in) 
   {
      readOutputFileName(in);
   }
   
   /*
   * OneLoopAnalyzer setup
   */
   template <int D>
   void OneLoopAnalyzer<D>::setup() 
   {
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);
   
      //Allocate variables
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      int meshSize = system().domain().mesh().size();
      
      if (!isInitialized_){
         fluctuation_.allocate(dimensions);
         
         // Allocate intraCorrelations_ S_AA, S_BB, SAB
         intraCorrelations_.allocate(3);
         for (int i = 0; i < 3; ++i) {
            intraCorrelations_[i].allocate(meshSize);
         }
         
         isInitialized_ = true;
      }
      
      computeIntra();
      computeFluctuation();
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }
   
   template <int D>
   void OneLoopAnalyzer<D>::computeIntra()
   {
      const double vMonomer = system().mixture().vMonomer();
      Mixture<D> const & mixture = system().mixture();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      int meshSize = system().domain().mesh().size();
      
      double phi, cPolymer, polymerLength;
      double length, lengthA, lengthB, ksq;
      double kuhn, kuhnA, kuhnB, eA, eB;
      int monomerId, monomerIdA, monomerIdB;
      int rank;
      
      IntVec<D> G, Gmin;
      DArray<double> Gsq;
      Gsq.allocate(meshSize);
      MeshIterator<D> iter;
      iter.setDimensions(dimensions);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, dimensions, unitCell);
         Gsq[iter.rank()] = unitCell.ksq(Gmin);
      }
      
      for (int i = 0; i < 3; ++i){
         for (int j = 0; j < meshSize; ++j){
            intraCorrelations_[i][j] = 0.0;
         }
      }
      
      PolymerSpecies const & polymer = system().mixture().polymerSpecies(0);
      phi = polymer.phi();
      polymerLength = 0.0;
      int nBlock = polymer.nBlock();
      for (int j = 0; j < nBlock; j++) {
         if (PolymerModel::isThread()) {
            length = polymer.edge(j).length();
         } else {
            length = (double) polymer.edge(j).nBead();
         }
         polymerLength += length;
      }
      
      // Compute cPolymer (polymer number concentration)
      cPolymer = phi/(polymerLength*vMonomer);
      
      // Compute diagonal SAA, SBB
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         ksq = Gsq[rank];
         
         // SAA, SBB
         for (int j = 0; j < nBlock; j++) {
            monomerId = polymer.edge(j).monomerId();
            kuhn = mixture.monomer(monomerId).kuhn();
            
            if (PolymerModel::isThread()) {
               length = polymer.edge(j).length();
               intraCorrelations_[j][rank] = cPolymer * Debye::dt(ksq, length, kuhn);
            } else {
               length = (double) polymer.edge(j).nBead();
               intraCorrelations_[j][rank] = cPolymer * Debye::db(ksq, length, kuhn);
            }
         }

         // SAB
         monomerIdA = polymer.edge(0).monomerId();
         kuhnA = mixture.monomer(monomerIdA).kuhn();
         monomerIdB = polymer.edge(1).monomerId();
         kuhnB = mixture.monomer(monomerIdB).kuhn();
         if (PolymerModel::isThread()) {
            lengthA = polymer.edge(0).length();
            lengthB = polymer.edge(1).length();
            eA = Debye::et(ksq, lengthA, kuhnA);
            eB = Debye::et(ksq, lengthB, kuhnB);
            intraCorrelations_[2][rank] = cPolymer * eA * eB;
         } else {
            lengthA = (double) polymer.edge(0).nBead();
            lengthB = (double) polymer.edge(1).nBead();
            eA = Debye::eb(ksq, lengthA, kuhnA);
            eB = Debye::eb(ksq, lengthB, kuhnB);
            intraCorrelations_[2][rank] = cPolymer * eA * eB;
         }
      }
   }
   
   /**
   * Compute fluctuation function
   */ 
   template <int D>
   void OneLoopAnalyzer<D>::computeFluctuation()
   {
      const double vMonomer = system().mixture().vMonomer();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      MeshIterator<D> iter;
      iter.setDimensions(dimensions);
      double Spp; double Smm; double Spm;
      double Saa; double Sbb; double Sab;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         Saa = vMonomer * intraCorrelations_[0][iter.rank()];
         Sbb = vMonomer * intraCorrelations_[1][iter.rank()];
         Sab = vMonomer * intraCorrelations_[2][iter.rank()];
         Spp = Saa + Sbb + 2.0* Sab;
         Smm = Saa + Sbb - 2.0* Sab;
         Spm = Saa - Sbb;
         fluctuation_[iter.rank()] = 4.0 * Spp / (Spp * Smm - Spm * Spm);
      }
      
   }
   
   /*
   * Output final results to output file.
   */
   template <int D>  
   void OneLoopAnalyzer<D>::output() 
   {
      double chi= system().interaction().chi(0,1);
      double freeEnergy = 0.0;
      int meshSize = system().domain().mesh().size();
      
      for (int i = 0; i < meshSize; i++){
         freeEnergy += 1.0/2.0 * std::log( 1.0 - 2.0* chi/fluctuation_[i]);
      }
      
      // Output structure factors to one file
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      
      outputFile_ << "One loop contribution is " <<freeEnergy <<" k_bT";
      outputFile_<< std::endl;
      outputFile_.close();
      Log::file() << "One loop contribution is " <<freeEnergy <<" k_bT";
   }

}
}
#endif 

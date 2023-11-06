#ifndef PSPC_COMPRESSOR_H
#define PSPC_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>   
#include <prdc/cpu/RField.h> 
#include <prdc/crystal/shiftToMinimum.h>   
#include <pscf/mesh/MeshIterator.h>     
#include <pspc/solvers/Polymer.h>

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Pspc_Compressor_Module
   */
   template <int D>
   class Compressor : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Compressor();

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Compressor(System<D>& system);

      /**
      * Destructor.
      */
      ~Compressor();

      /**
      * Iterate Langrange multiplier field.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int compress() = 0;
      
      /**
      * Count how many times MDE have been computed
      */
      virtual int counterMDE() = 0;
      
      /**
      * Return Anderson mixing  subspace step reduce error percentage
      */
      virtual double subspacePercent() = 0;
      
      /**
      * Return Anderson mixing correction step reduce error percentage
      */
      virtual double correctionPercent() = 0;
      
      /**
      * Log output timing results 
      */
      virtual void outputTimers(std::ostream& out) = 0;
      
      /**
      * Clear timers 
      */
      virtual void clearTimers() = 0;
      
      /**
      * Compute Debye function
      */
      double computeDebye(double x);
      
      /**
      * Compute intramolecular correlation at specific sqSquare
      */
      double computeIntraCorrelation(double qSquare);
      
      /**
      * Compute intramolecular correlation  
      */
      RField<D> computeIntraCorrelation();
      
      /**
      * Return const reference to parent system.
      */
      System<D> const & system() const
      {  return *sysPtr_;}

   protected:

      /**
      * Return reference to parent system.
      */
      System<D>& system() 
      {  return *sysPtr_;}

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Inline member functions

   // Default constructor
   template <int D>
   inline Compressor<D>::Compressor()
    : sysPtr_(&system)
   {  setClassName("Compressor"); }

   // Constructor
   template <int D>
   Compressor<D>::Compressor(System<D>& system)
    : sysPtr_(&system)
   {  setClassName("Compressor"); }

   // Destructor
   template <int D>
   Compressor<D>::~Compressor()
   {}
   
   template<int D>
   double Compressor<D>::computeDebye(double x)
   {
      if (x == 0){
         return 1.0;
      } else {
         return 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
   }
   
   template<int D>
   double Compressor<D>::computeIntraCorrelation(double qSquare)
   {
      const int np = system().mixture().nPolymer();
      const double vMonomer = system().mixture().vMonomer();
      // Overall intramolecular correlation
      double omega = 0;
      int monomerId; int nBlock; 
      double kuhn; double length; double g; double rg2; 
      Polymer<D> const * polymerPtr;
      
      for (int i = 0; i < np; i++){
         polymerPtr = &system().mixture().polymer(i);
         nBlock = polymerPtr->nBlock();
         double totalLength = 0;
         for (int j = 0; j < nBlock; j++) {
            totalLength += polymerPtr->block(j).length();
         }
         for (int j = 0; j < nBlock; j++) {
            monomerId = polymerPtr-> block(j).monomerId();
            kuhn = system().mixture().monomer(monomerId).kuhn();
            // Get the length (number of monomers) in this block.
            length = polymerPtr-> block(j).length();
            rg2 = length * kuhn* kuhn /6.0;
            g = computeDebye(qSquare*rg2);
            omega += length/totalLength * length * g/ vMonomer;
         }
      }
      return omega;
   }
   
   template<int D>
   RField<D> Compressor<D>::computeIntraCorrelation()
   {
      RField<D> intraCorrelation;
      IntVec<D> const & dimensions = system().mesh().dimensions();
      IntVec<D> kMeshDimensions;
      MeshIterator<D> iter;
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions[i] = dimensions[i];
         } else {
            kMeshDimensions[i] = dimensions[i]/2 + 1;
         }
      }
      iter.setDimensions(kMeshDimensions);
      intraCorrelation.allocate(kMeshDimensions);
      IntVec<D> G, Gmin;
      double Gsq;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
         Gsq = system().unitCell().ksq(Gmin);
         intraCorrelation[iter.rank()] = computeIntraCorrelation(Gsq);
      }
      return intraCorrelation;
   }
   
} // namespace Pspc
} // namespace Pscf
#endif

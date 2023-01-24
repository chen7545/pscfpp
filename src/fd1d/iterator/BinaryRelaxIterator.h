#ifndef FD1D_BINARY_RELAX_ITERATOR_H
#define FD1D_BINARY_RELAX_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <fd1d/solvers/Mixture.h>
#include <pscf/math/LuSolver.h>
#include <util/containers/Array.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>


namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /**
   * Relaxation iterator for SCF equations for two-monomer system.
   *
   * This class implements the simple Picard-type relaxation iterator 
   * for systems with two monomer types that was introduced by Drolet and 
   * Fredrickson (PRL, 1999).
   * 
   * Reference:
   * F. Drolet & G.H. Fredrickson, Phys. Rev. Lett. vol. 83, 4317 (1999).
   *
   * \ingroup Fd1d_Iterator_Module
   */
   class BinaryRelaxIterator : public Iterator
   {
    
   public:
    
      /**
      * Monomer chemical potential field.
      */
      typedef Mixture::WField WField;
    
      /**
      * Monomer concentration / volume fraction field.
      */
      typedef Mixture::CField CField;
   
      /**
      * Constructor.
      *
      * \param system parent System object.
      */
      BinaryRelaxIterator(System& system);
    
      /**
      * Destructor.
      */
      virtual ~BinaryRelaxIterator();
    
      /**
      * Read all parameters and initialize.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);
    
      /**
      * Iterate self-consistent field equations to solution.
      *
      * \param isContinuation True if part of sweep, and not first step.
      * \return error code: 0 for success, 1 for failure.
      */
      int solve(bool isContinuation = false);
   
   private:
    
    
      /// Perturbed chemical potential fields (work space).
      DArray<WField> wFieldsNew_;
    
      /// Perturbed monomer concentration fields (work space).
      DArray<WField> cFieldsNew_;
      
      
      /// Concentrations at one point (work space).
      DArray<double> cArray_;
    
      /// Chemical potentials at one point (work space).
      DArray<double> wArray_;
    
      /// Residual vector. size = nr = (# monomers)x(# grid points).
      DArray<double> residual_;
    
      /// Perturbed residual. size = nr.
      DArray<double> residualNew_;
      
      /// Perturbed DW
      DArray<WField> dW_;
      
      DArray<WField> dWNew_;
      
      /// Norm of change field
      double dWNorm_;
      
      double dWNormNew_;
      
      /// Error tolerance.
      double epsilon_;
      
      /// Mixing parameter for Wplus 
      double lambdaPlus_;
      
      /// Mixing parameter for Wminus
      double lambdaMinus_;
      
      /// Max number of iterations
      int maxIter_;
    
      /// Have arrays been allocated?
      bool isAllocated_;
    
      /// Is the ensemble canonical for all species ?
      bool isCanonical_;
      
   
      /**
      * Allocate required memory.
      */
      void allocate();
    
      /**
      * Compute residuals and increments of chemical potential fields.
      *
      * \param wOld  array of old chemical potential fields (input)
      * \param cFields  monomer concentration fields (input)
      * \param dW  change of w fields (output)
      * \param dWNorm scalar residual (output)
      */
      void computeDW(Array<WField> const & wOld, 
                     Array<CField> const & cFields,
                     Array<WField> & dW,
                     double & dWNorm);
     
      /**
      * Update the chemical potential fields
      *
      * \param wold array of old chemical potential fields (input)
      * \param dW_ array of increment of chemical potential fields (input)
      * \param wNew array of new chemical potential fields (output)
      */
      void updateWFields(Array<WField> const & wOld,
                         Array<WField> const & dW_,
                         Array<WField> & wNew);
                         
      void computeResidual(Array<WField> const & wFields, 
                           Array<WField> const & cFields, 
                           Array<double>& residual);
                           
      double residualNorm(Array<double> const & residual) const;
                           
                         
   };

} // namespace Fd1d
} // namespace Pscf
#endif

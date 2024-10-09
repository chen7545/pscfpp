#ifndef RPG_AM_COMPRESSOR_H
#define RPG_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"
#include <prdc/cuda/Field.h>
#include <prdc/cuda/RField.h>         

#include <pscf/iterator/AmIteratorTmpl.h>     

#include <util/containers/DArray.h>                 

namespace Pscf {
namespace Rpg
{

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Rpg implementation of the Anderson Mixing compressor.
   *
   * \ingroup Rpg_Simulate_Compressor_Module
   */
   template <int D>
   class AmCompressor : public AmIteratorTmpl<Compressor<D>, Field<cudaReal> >
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this compressor.
      */
      AmCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~AmCompressor();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at the 
      * beginning of iteration
      *
      * \param isContinuation true iff continuation within a sweep
      */ 
      void setup(bool isContinuation);      
      
      
      /**
      * Compress to obtain partial saddle point w+
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress();    
      
      /**
      * Return how many times MDE has been solved.
      */
      int mdeCounter(); 
      
      /**
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out);
      
      /**
      * Clear all timers (reset accumulated time to zero).
      */
      void clearTimers();
      
      /**
      * Compute each step of error
      */
      double computeInCompressError();
      
      /**
      * Get the ratio of error reduction by AM step 1
      */
      std::vector<double> stepOneRatioVector();
      
      /**
      * Get the predicted of error reduction by AM step 1
      */
      std::vector<double> predictRatioVector();
      
      /**
      * Get the ratio of error reduction by AM step 2
      */
      std::vector<double> stepTwoRatioVector();
      
      // Inherited public member functions
      using AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::setClassName;

   protected:
  
      // Inherited protected members 
      using ParamComposite::readOptional;
      using Compressor<D>::system;
      using Compressor<D>::mdeCounter_;
      using Compressor<D>::totalItr_;

   private:
      
      std::vector<double> stepOneRatioVector_; 
      std::vector<double> predictRatioVector_; 
      std::vector<double> stepTwoRatioVector_; 
      
      /**
      * Type of error criterion used to test convergence 
      */ 
      std::string errorType_;
      
      /**
      * Incompressibility constraint error.
      */ 
      RField<D> error_;
      
      /**
      * Count how many times MDE has been solved.
      */
      int counter_;
      
      /**
      * Current values of the fields
      */
      DArray< RField<D> > w0_;  

      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
      
      /**
      * Template w Field used in update function
      */
      DArray< RField<D> > wFieldTmp_;
      
      /**
      * New Basis variable used in updateBasis function 
      */
      Field<cudaReal> newBasis_;

      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(Field<cudaReal>& a, Field<cudaReal> const & b);

      /**
      * Compute the inner product of two vectors
      */
      double dotProduct(Field<cudaReal> const & a, Field<cudaReal> const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(Field<cudaReal> const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<Field<cudaReal> > & basis, 
                       RingBuffer<Field<cudaReal> > const & hists);

      /**
      * Add linear combination of basis vectors to trial field.
      * 
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(Field<cudaReal>& trial, 
                        RingBuffer<Field<cudaReal> > const & basis, 
                        DArray<double> coeffs, 
                        int nHist);

      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(Field<cudaReal>& fieldTrial, 
                             Field<cudaReal> const & resTrial, 
                             double lambda);

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess();
     
      /** 
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements();

      /**
      * Gets the current field vector from the system.
      * 
      * \param curr current field vector
      */ 
      void getCurrent(Field<cudaReal>& curr);

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate();

      /**
      * Compute the residual vector.
      *
      * \param resid current residual vector value
      */
      void getResidual(Field<cudaReal>& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(Field<cudaReal>& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();
      

   };
   
   // Inline functions

   // Get the ratio of error reduction by AM step 1
   template <int D>
   inline std::vector<double> AmCompressor<D>::stepOneRatioVector()
   { return stepOneRatioVector_; }
   
   // Get the ratio of error reduction by AM step 2
   template <int D>
   inline std::vector<double> AmCompressor<D>::stepTwoRatioVector()
   { return stepTwoRatioVector_; }
   
   // Get the ratio of error reduction by AM step 2
   template <int D>
   inline std::vector<double> AmCompressor<D>::predictRatioVector()
   { return predictRatioVector_; }
   
} // namespace Rpg
} // namespace Pscf
#endif

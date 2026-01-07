#ifndef PSCF_CPU_REDUCE_H
#define PSCF_CPU_REDUCE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Functions that perform array reductions on a CPU.
   *
   * A reduction is any operation that involves reducing all of 
   * the elements of an array or set of arrays to a single scalar.  
   * Examples include taking the sum or finding the maximum of all 
   * array elements, or taking an inner product of two arrays.
   *
   * \defgroup Pscf_Cpu_Reduce_Module Reduce (CPU)
   * \ingroup Pscf_Cpu_Module
   */

   namespace Reduce {

      /*
      *  This file and Reduce.cpp contain declarations and definitions
      *  for reduction functions that operate real arrays.  Corresponding 
      *  declarations and definitions for reductions of complex-valued
      *  arrays are in files ReduceCx.h and ReduceCx.cpp, respectively. 
      */

      // Summation

      /**
      * Compute sum of array elements (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double sum(Array<double> const & in);

      /**
      * Compute sum of of squares of array elements (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double sumSq(Array<double> const & in);

      // Summation of products

      /**
      * Compute Euclidean inner product of two real arrays .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param a  first input array
      * \param b  second input array
      */
      double innerProduct(Array<double> const & a,
                          Array<double> const & b);

      // Maxima and minima

      /**
      * Get maximum of array elements .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double max(Array<double> const & in);

      /**
      * Get maximum absolute magnitude of array elements .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double maxAbs(Array<double> const & in);

      /**
      * Get minimum of array elements .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double min(Array<double> const & in);

      /**
      * Get minimum absolute magnitude of array elements .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double minAbs(Array<double> const & in);

   } // namespace Reduce
} // namespace Pscf
#endif

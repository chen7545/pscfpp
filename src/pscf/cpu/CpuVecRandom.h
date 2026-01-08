#ifndef PSCF_CPU_VEC_RANDOM_H
#define PSCF_CPU_VEC_RANDOM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declarations
namespace Util {
   class Random;
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Random number generator for arrays of random numbers on CPU.
   *
   * This class generates arrays of random numbers on a CPU,
   * using an interface analogous to that of CudaRandom to allow
   * creation of template code that can use either type of object.
   * It uses an associated Util::Random object to generate random
   * numbers.
   *
   * \ingroup Pscf_Cpu_Module
   */
   class CpuVecRandom
   {

   public:

      /**
      * Default constructor.
      */
      CpuVecRandom();

      /**
      * Constructor - creates association
      *
      * \param random  associated scalar random number generator
      */
      CpuVecRandom(Util::Random& random);

      /**
      * Destructor.
      */
      virtual ~CpuVecRandom();

      /**
      * Create an association with a Util::Random object.
      *
      * \param random  associated scalar random number generator
      */
      void associate(Util::Random& random);

      /**
      * Populate array on device with random doubles in (0, 1], uniform dist.
      *
      * \param data  array to populate
      */
      void uniform(Array<double>& data);

      /**
      * Populate array on device with normal-distributed random doubles.
      *
      * Note: the input array must have an even number of elements. This is a
      * requirement imposed by cuRAND, the random number generator software
      * used by CpuVecRandom.
      *
      * \param data  array to populate
      * \param stddev  standard deviation (input)
      * \param mean  mean value (input, default = 0.0)
      */
      void normal(Array<double>& data, double stddev, double mean = 0.0);

   private:

      /// Pointer to associated scalar random number generator.
      Util::Random* randomPtr_;

   };

   #if 0
   /*
   * Returns value of random seed (private member variable idum)
   */
   inline long CpuVecRandom::seed()
   {  return randomPtr_->seed(); }
   #endif

}
#endif

#include "CpuRandom.h"
#include <util/random/Random.h>
#include <util/containers/Array.h>
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   CpuRandom::CpuRandom()
    : randomPtr_(nullptr)
   {}

   /*
   * Associating constructor.
   */
   CpuRandom::CpuRandom(Util::Random& random)
    : randomPtr_(&random)
   {}

   /*
   * Destructor.
   */
   CpuRandom::~CpuRandom()
   {}

   /*
   * Create an association after construction
   */
   void CpuRandom::associate(Util::Random& random)
   {  randomPtr_ = &random; }

   /*
   * Populate array on device with random doubles in (0, 1], uniform dist.
   */
   void CpuRandom::uniform(Array<double>& data)
   {
      UTIL_CHECK(randomPtr_);
      UTIL_CHECK(data.capacity() > 0);
      const int n = data.capacity();
      for (int i = 0; i < n; ++i) {
         data[i] = randomPtr_->uniform();
      } 
   }

   /*
   * Populate array with normal-distributed random doubles.
   */
   void CpuRandom::normal(Array<double>& data, 
                          double stddev, double mean)
   {
      UTIL_CHECK(randomPtr_);
      UTIL_CHECK(data.capacity() > 0);
      const int n = data.capacity();
      for (int i = 0; i < n; ++i) {
         data[i] = mean + stddev * randomPtr_->gaussian();
      }
   }

} // namespace Pscf

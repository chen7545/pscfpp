/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceMemory.h"
#include "cudaErrorCheck.h"
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   DeviceMemory::DeviceMemory()
    : dataPtr_(nullptr),
      capacity_(0)
   {}

   /*
   * Allocating constructor.
   */
   DeviceMemory::DeviceMemory(int capacity)
    : dataPtr_(nullptr),
      capacity_(0)
   {  allocate(capacity); }

   /*
   * Destructor.
   */
   DeviceMemory::~DeviceMemory()
   {
      if (isAllocated()) {
         cudaFree(dataPtr_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   void DeviceMemory::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate DeviceMemory");
      }
      cudaErrorCheck( cudaMalloc( &dataPtr_, capacity ) );
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this array is not allocated.
   */
   void DeviceMemory::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Attempt to deallocate unallocated DeviceMemory");
      }
      cudaErrorCheck( cudaFree(dataPtr_) );
      dataPtr_ = nullptr;
      capacity_ = 0;
   }

   /*
   * Reallocate if necessary to increase capacity.
   */
   void DeviceMemory::resize(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to resize DeviceMemory with capacity <= 0");
      }
      if (capacity > capacity_) {
         if (isAllocated()) {
            deallocate();
         }
         allocate(capacity);
      }
   }

   /*
   * Get a pointer to the underlying C array.
   */
   void* DeviceMemory::cArray() const
   {
      UTIL_CHECK(dataPtr_);
      return dataPtr_;
   }

} // namespace Pscf

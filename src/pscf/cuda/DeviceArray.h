#ifndef PSCF_DEVICE_ARRAY_H
#define PSCF_DEVICE_ARRAY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/ReferenceCounter.h>
#include <util/misc/CountedReference.h>

namespace Pscf {

   // Forward declarations
   template <typename Data> class HostDArray;
   class DeviceMemory;

   using namespace Util;

   /**
   * Dynamic array on the GPU device with aligned data.
   *
   * This class wraps an aligned C array with elements of type Data that
   * is allocated in GPU device global memory.  All member functions may
   * be called from the CPU host, but this class does not offer access
   * to individual elements via the subscript operator, operator[].
   *
   * A DeviceArray can have one of two different relationships with its
   * underlying data. In case 1, the data are owned by this object, so
   * the allocation and destruction of the C array are performed by this
   * object. In case 2, this object is instead "associated" with a slice
   * of an array owned by a different DeviceArray. If this is the case,
   * this object is not responsible for allocation or destruction of the
   * underlying C array, and merely acts as a reference to the slice of
   * the other array to which it is associated.
   *
   * \ingroup Pscf_Cuda_Module
   */
   template <typename Data>
   class DeviceArray
   {

   public:

      /**
      * Data type of each element.
      */
      typedef Data ValueType;

      /**
      * Default constructor.
      */
      DeviceArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      *
      * \param capacity number of elements to allocate
      */
      DeviceArray(int capacity);

      /**
      * Copy constructor.
      *
      * \param other DeviceArray<Data> to be copied (input)
      */
      DeviceArray(DeviceArray<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DeviceArray();

      /**
      * Allocate the underlying C array on the device.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the array is not allocated.
      * \throw Exception if this object does not own the array.
      */
      void deallocate();

      /**
      * Associate this object with a slice of a different DeviceArray.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param arr  parent array that owns the data
      * \param beginId  index in the parent array at which this array starts
      * \param capacity  number of elements associated with this container
      */
      void associate(DeviceArray<Data>& arr, int beginId, int capacity);

      /**
      * Associate this object with a DeviceMemory container.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param arr  DeviceMemory container that owns the data
      * \param capacity  number of elements of type Data
      */
      void associate(DeviceMemory& arr, int capacity);

      /**
      * Dissociate this object from an associated external memory block.
      *
      * \throw Exception if this object is not associated external memory.
      */
      void dissociate();

      /**
      * Associate a reference with reference counter for this container.
      *
      * This function should normally be not called by users of a
      * DeviceArray container. It is usually called within the associate
      * member function of a container that is referring to data owned
      * by this array.
      *
      * \param reference  reference to be included by reference counter
      */
      void addReference(CountedReference& reference);

      /**
      * Assignment operator, assign from another DeviceArray<Data> array.
      *
      * Performs a deep copy, by copying values of all elements from
      * device memory to device memory.
      *
      * This function will allocate memory if this (LHS) array is not
      * allocated.  If this array is arleady allocated, it must have the
      * same capacity as the other (RHS) DeviceArray<Data>.
      *
      * \param other DeviceArray<Data> on rhs of assignent (input)
      */
      virtual
      DeviceArray<Data>& operator = (const DeviceArray<Data>& other);

      /**
      * Assignment operator, assignment from HostDArray<Data> host array.
      *
      * Performs a deep copy from a RHS HostDArray<Data> host array to
      * this LHS DeviceArray<Data> device array, by copying underlying
      * C array from host memory to device memory.
      *
      * This function will allocate memory if this (LHS)
      * DeviceArray<Data> is not allocated.  If this is allocated, it
      * must have the same dimensions as the RHS HostDArray<Data>.
      *
      * \param other HostDArray<Data> on RHS of assignent (input)
      */
      virtual
      DeviceArray<Data>& operator = (const HostDArray<Data>& other);

      /**
      * Return array capacity.
      *
      * \return number of elements in array
      */
      int capacity() const;

      /**
      * Return true if the array has allocated data, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Does this container have an array that it owns?
      */
      bool isOwner() const;

      /**
      * Is this container associated with data it does not own?
      */
      bool isAssociated() const;

      /**
      * Return pointer to underlying C array.
      */
      Data* cArray();

      /**
      * Return const pointer to underlying C array.
      */
      Data const * cArray() const;

   protected:

      /// Pointer to a C array of Data elements on the GPU device.
      Data* dataPtr_;

      /// Allocated size (capacity) of the array.
      int capacity_;

      /// Counter for any arrays that use data owned by this.
      ReferenceCounter refCounter_;

      /// Reference to another array that owns memory used by this one.
      CountedReference ref_;

   };

   // Inline functions

   /*
   * Return array capacity.
   */
   template <typename Data>
   inline int DeviceArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return true if the array has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool DeviceArray<Data>::isAllocated() const
   {  return (bool) dataPtr_; }

   /*
   * Does this object own data?
   */
   template <typename Data>
   inline bool DeviceArray<Data>::isOwner() const
   {  return ((bool) dataPtr_ && !ref_.isAssociated()); }

   /*
   * Is this object associated with data it does not own?
   */
   template <typename Data>
   inline bool DeviceArray<Data>::isAssociated() const
   {  return (bool) dataPtr_ && ref_.isAssociated(); }

}

#include "DeviceMemory.h"
#include "HostDArray.h"
#include "cudaErrorCheck.h"
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Default constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray()
    : dataPtr_(nullptr),
      capacity_(0)
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(int capacity)
    : dataPtr_(nullptr),
      capacity_(0)
   {  allocate(capacity); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(const DeviceArray<Data>& other)
    : dataPtr_(nullptr),
      capacity_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other array must be allocated.");
      }

      allocate(other.capacity_);
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyDeviceToDevice) );
   }

   /*
   * Destructor.
   */
   template <typename Data>
   DeviceArray<Data>::~DeviceArray()
   {
      if (isOwner()) {
         if (refCounter_.hasRefs()) {
            std::cout
                << "Error: Destroying DeviceArray that is referred"
                << "to by at least one CountedReference"
                << std::endl;
         }
         cudaFree(dataPtr_);
      }
      if (ref_.isAssociated()) {
         ref_.dissociate();
      }
   }

   /*
   * Allocate the underlying C array.
   */
   template <typename Data>
   void DeviceArray<Data>::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate an array");
      }
      if (ref_.isAssociated()) {
         UTIL_THROW("Attempt to allocate array that references other data");
      }

      cudaErrorCheck(cudaMalloc((void**) &dataPtr_, capacity * sizeof(Data)));
      capacity_ = capacity;

   }

   /*
   * Deallocate the underlying C array, if possible.
   */
   template <typename Data>
   void DeviceArray<Data>::deallocate()
   {
      UTIL_CHECK(dataPtr_);
      UTIL_CHECK(!ref_.isAssociated());
      UTIL_CHECK(!refCounter_.hasRefs());
      cudaErrorCheck( cudaFree(dataPtr_) );
      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
   }

   /*
   * Associate this object with a slice of a different DeviceArray.
   */
   template <typename Data>
   void DeviceArray<Data>::associate(DeviceArray<Data>& arr, int beginId,
                                     int capacity)
   {
      UTIL_CHECK(arr.isAllocated());
      UTIL_CHECK(arr.isOwner());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(beginId + capacity <= arr.capacity());
      UTIL_CHECK(!dataPtr_);
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and capacity
      dataPtr_ = arr.cArray() + beginId;
      capacity_ = capacity;

      // Associate ref_ member with reference counter of data owner array
      arr.addReference(ref_);
   }

   /*
   * Associate this object with a slice of a different DeviceArray.
   */
   template <typename Data>
   void DeviceArray<Data>::associate(DeviceMemory& arr, int capacity)
   {
      UTIL_CHECK(arr.isAllocated());
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(capacity * sizeof(Data) <= arr.capacity());
      UTIL_CHECK(!dataPtr_);
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and capacity
      dataPtr_ = arr.cArray();
      capacity_ = capacity;

      // Associate ref_ member with reference counter of data owner array
      arr.addReference(ref_);
   }

   /*
   * Dissociate this object from the associated array.
   *
   * Throw Exception if this object is not associated with another array.
   */
   template <typename Data>
   void DeviceArray<Data>::dissociate()
   {
      UTIL_CHECK(dataPtr_);
      UTIL_CHECK(ref_.isAssociated());

      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
      ref_.dissociate();
   }

   /*
   * Associate a CountedReference with the ReferenceCounter.
   */
   template <typename Data>
   void DeviceArray<Data>::addReference(CountedReference& ref)
   {  ref.associate(refCounter_); }

   /*
   * Assignment from another DeviceArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>&
   DeviceArray<Data>::operator = (const DeviceArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - RHS DeviceArray<Data> must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("Other DeviceArray<Data> must be allocated.");
      }

      // If this is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Require equal capacity values
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy elements
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyDeviceToDevice) );

      return *this;
   }

   /*
   * Assignment of LHS DeviceArray<Data> from RHS HostDArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>&
   DeviceArray<Data>::operator = (const HostDArray<Data>& other)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("RHS HostDArray<Data> must be allocated.");
      }

      // Allocate this if necessary
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Require equal capacity values
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy elements
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyHostToDevice) );

      return *this;
   }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   Data* DeviceArray<Data>::cArray()
   {
      UTIL_CHECK(dataPtr_);
      return dataPtr_;
   }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   const Data* DeviceArray<Data>::cArray() const
   {
      UTIL_CHECK(dataPtr_);
      return dataPtr_;
   }

} // namespace Pscf
#endif

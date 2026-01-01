#ifndef PSCF_DEVICE_MEMORY_H
#define PSCF_DEVICE_MEMORY_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   /**
   * Block of bare memory allocated on the GPU device.
   *
   * This class wraps a void* C array that is allocated in GPU device 
   * global memory, and also contains the number of bytes allocated.
   *
   * \ingroup Pscf_Cuda_Module
   */
   class DeviceMemory
   {

   public:

      /**
      * Default constructor.
      */
      DeviceMemory();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity size in bytes to allocate 
      */
      DeviceMemory(int capacity);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DeviceMemory();

      /**
      * Allocate the underlying C array on the device.
      *
      * \throw Exception if capacity is not positve
      * \throw Exception if the array is already allocated.
      *
      * \param capacity  number of byes to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the array is not allocated.
      */
      void deallocate();

      /**
      * Re-allocate if necessary to increase capacity.
      *
      * If the capacity parameter is greater than current capacity, allocate
      * or re-allocate to obtain a block of desired capacity. Otherwise, do
      * nothing. 
      *
      * When a previously allocated array is resized to increase capacity,
      * the old array is deallocated and all data is lost. 
      *
      * \throw Exception if capacity is not positve
      *
      * \param capacity  number of byes requird
      */
      void resize(int capacity);

      /**
      * Return pointer to underlying C array.
      */
      void* cArray() const;

      /**
      * Return allocated capacity.
      *
      * \return Number of elements allocated in array
      */
      int capacity() const;

      /**
      * Return true if the array has been allocated, false otherwise.
      */
      bool isAllocated() const;

   protected:

      /// Pointer to a memory block on the GPU device.
      void* dataPtr_;

      /// Allocated size (capacity) of the array in bytes.
      int capacity_;

   };

   /*
   * Return true if the array has been allocated, false otherwise.
   */
   inline
   bool DeviceMemory::isAllocated() const
   {  return (bool)dataPtr_; }

   /*
   * Return allocated capacity.
   */
   inline 
   int DeviceMemory::capacity() const
   {  return capacity_; }

} // namespace Pscf
#endif

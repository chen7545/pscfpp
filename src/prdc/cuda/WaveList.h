#ifndef PRDC_CUDA_WAVE_LIST_H
#define PRDC_CUDA_WAVE_LIST_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RField.h>        // member
#include <pscf/cuda/DeviceArray.h>   // member
#include <pscf/cuda/HostDArray.h>    // member
#include <pscf/math/IntVec.h>        // member
#include <util/containers/DArray.h>  // member

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Class to calculate and store properties of wavevectors.
   *
   * In particular, minimum images, square norms of wavevectors (kSq), and
   * derivatives of the square norms of wavevectors with respect to the
   * lattice parameters (dKSq) are calculated and stored by this class.
   *
   * Any time the lattice parameters change the clearUnitCellData() method
   * should be called, which will effectively reset the WaveList object so
   * that the wavevector properties will need to be recalculated before
   * being used.
   *
   * This object calculates these wavevector properties for a mesh of grid
   * points in k-space. If a calculation only requires real-valued fields,
   * PSCF uses a reduced-size k-space mesh, as output by cuFFT. However, a
   * full-sized k-space mesh (the same size as the real-space mesh) is
   * necessary when dealing with complex-valued fields. The k-space mesh
   * used by a WaveList object is determined by the parameter isRealField,
   * which is assigned in the constructor and cannot later be changed.
   */
   template <int D>
   class WaveList
   {
   public:

      /**
      * Constructor.
      *
      * \param isRealField  Will this WaveList be used for real-valued fields?
      */
      WaveList(bool isRealField = true);

      /**
      * Destructor
      */
      ~WaveList();

      /**
      * Allocate memory and set association with a Mesh and UnitCell object.
      *
      * \param m  spatial discretization mesh (input)
      * \param c  crystallographic unit cell (input)
      */
      void allocate(Mesh<D> const & m, UnitCell<D> const & c);

      /**
      * Clear all internal data that depends on lattice parameters.
      *
      * Sets hasKSq_ and hasdKSq_ to false, and sets hasMinImages_ to
      * false if the crystall lattice type has variable angles.
      */
      void clearUnitCellData();

      /**
      * Compute minimum images of wavevectors. (Also calculates kSq.)
      *
      * The minimum images may change if a lattice angle in the unit cell
      * is changed, so this method should be called whenever such changes
      * occur.
      *
      * In the process of computing the minimum images, the square norm
      * |k|^2 for all wavevectors is also calculated and stored, so it
      * is not necessary to call computeKSq after calling this method.
      * computeKSq is provided to allow calculation of kSq without
      * recalculating minimum images.
      */
      void computeMinimumImages();

      /**
      * Compute square norm |k|^2 for all wavevectors.
      *
      * This function uses existing mininum images if they are valid, or
      * recomputes them if necessary.
      */
      void computeKSq();

      /**
      * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
      */
      void computedKSq();

      /**
      * Get the array of minimum images on the device by reference.
      *
      * The array has size kSize * D, where kSize is the number of grid
      * points in the FFT k-space mesh. The array is unwrapped into a
      * linear array in an index-by-index manner, in which the first kSize
      * elements of the array contain the first index of each minimum
      * image, and so on. If isRealField is true, kSize is smaller than
      * the size of the real-space mesh. Otherwise, it is equal.
      */
      DeviceArray<int> const & minImages_d() const;

      /**
      * Get minimum images as IntVec<D> objects on the host.
      *
      * The array has size kSize, and each element is an IntVec<D>.
      * If isRealField is true, kSize is smaller than the size of the
      * real-space mesh. Otherwise, it is equal.
      */
      HostDArray< IntVec<D> > const & minImages_h() const;

      /**
      * Get the kSq array on the device by reference.
      *
      * This method returns an RField in which each element is the square
      * magnitude |k|^2 of a wavevector k in the k-space mesh used for the
      * DFT. If isRealField is true, this k-space mesh is smaller than the
      * real-space mesh. Otherwise, it is the same size.
      */
      RField<D> const & kSq() const;

      /**
      * Get derivatives of |k|^2 with respect to lattice parameter i.
      *
      * This method returns an RField in which each element is the
      * derivative of the square-wavevector with respect to unit cell
      * parameter i, multiplied by a prefactor. The prefactor is 2.0 for
      * waves that have an implicit inverse and 1.0 otherwise. The choice
      * of prefactor is designed to simplify use of the array to compute
      * stress.
      *
      * Each element corresponds to one wavevector k in the k-space mesh
      * used for the DFT. If isRealField is true, this k-space mesh is
      * smaller than the real-space mesh. Otherwise, it is the same size.
      * In the latter case, there are no implicit waves, so the prefactor
      * is always 1.0.
      *
      * \param i index of lattice parameter
      */
      RField<D> const & dKSq(int i) const;

      /**
      * Get the implicitInverse array by reference.
      *
      * This array is defined on a k-grid mesh, with a boolean value for
      * each gridpoint. The boolean represents whether the inverse of the
      * wave at the given gridpoint is an implicit wave. Implicit here is
      * used to mean any wave that is outside the bounds of the k-grid.
      *
      * This method will throw an error if isRealField == false, because
      * there are no implicit inverses in such a case.
      */
      DeviceArray<bool> const & implicitInverse() const;

      /**
      * Has memory been allocated for arrays?
      */
      bool isAllocated() const
      {  return isAllocated_; }

      /**
      * Have minimum images been computed?
      */
      bool hasMinImages() const
      {  return hasMinImages_; }

      /**
      * Has the kSq array been computed?
      */
      bool hasKSq() const
      {  return hasKSq_; }

      /**
      * Has the dKSq array been computed?
      */
      bool hasdKSq() const
      {  return hasdKSq_; }

      /**
      * Does this WaveList correspond to real-valued fields?
      */
      bool isRealField() const
      {  return isRealField_; }

   private:

      // Private member variables

      /**
      * Array containing minimum images for each wave, stored on device.
      *
      * The array has size kSize_ * D, where kSize_ is the number of grid
      * points in reciprocal space. The array is unwrapped into a linear
      * array in which the first kSize_ elements of the array contain
      * the the first coordinate for all minimum image, and so on. If
      * isRealField_ is true, kSize_ is smaller than the size of the
      * real-space mesh. Otherwise, kSize_ is equal to the size of the
      * real-space mesh.
      */
      DeviceArray<int> minImages_;

      /**
      * Array of IntVec<D> minimum images, stored on the host.
      *
      * Each element of minImageVecs_ contains all D coordinates of the
      * minimum image for a single wavevector, stored on the host as an
      * IntVec<D>. The array has capacity kSize_. If isRealField is true,
      * kSize_ is smaller than the size of the real-space mesh. Otherwise,
      * kSize_ is equal to the size of the real space mesh.
      */
      mutable
      HostDArray< IntVec<D> > minImages_h_;

      /**
      * Array containing values of kSq_, stored on the device.
      *
      * The mesh dimensions are those of the reciprocal space mesh.
      */
      RField<D> kSq_;

      /**
      * Array containing all values of dKSq_, stored on the device.
      *
      * The dimensions are kSize_ * nParam, where nParam is the number
      * of unit cell parameters.
      */
      DeviceArray<cudaReal> dKSq_;

      /**
      * Array of RFields, where each RField is a slice of the dKSq_ array.
      *
      * The number of elements is equal to nParam, the number of unit cell
      * parameters.  Element dKSqSlices_[i] is an  RField<D> element that
      * is associated with a slice of the larger dKSq_ device array, and
      * that contains derivatives of square wavevectors with respect to
      * unit cell parameter number i.
      *
      * The dKSqSlices_ container should appear after dKSq_ in the
      * declaration of class members in order to gurantee that the
      * elements of dkSqSlices_ will be destroyed before the dKSq_
      * container that owns the data.
      */
      DArray< RField<D> > dKSqSlices_;

      /**
      * Array indicating whether a given gridpoint has an implicit partner.
      *
      * This array is only allocated and used if isRealField_ is true,
      * in which case it has size kSize_.
      */
      DeviceArray<bool> implicitInverse_;

      /**
      * Dimensions of the mesh in reciprocal space.
      *
      * If isRealField_ is true, kMeshDimensions_ is equal to the vector
      * of dimensions of the reciprocal space grid used by a RFieldDft<D>
      * container to store the discrete Fourier transform of a real field.
      * One dimension of this mesh is approximately half the corresponding
      * dimension of the associated real space grid.  If isRealField_
      * is false, indicating application to a complex field, then
      * kMeshDimensions_ is equal to the vector of dimensions of the real
      * space grid.
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Number of grid points in reciprocal space.
      *
      * The integer kSize_ is the number of elements in the reciprocal
      * space grid, given by the product of elements of kMeshDimensions_.
      * If isRealField_, kSize_ is smaller than the size of the real
      * space mesh, by approximately a factor of 2 for large meshes.
      */
      int kSize_;

      /// Has memory been allocated for private member arrays?
      bool isAllocated_;

      /// Do valid minimum images exist (array minImages_) ?
      bool hasMinImages_;

      /// Have minimum image vectors been re-ordered in minImageVecs_ ?
      mutable
      bool hasMinImages_h_;

      /// Do valid values of kSq_ array exist ?
      bool hasKSq_;

      /// Do valid values of dKSq_ exist ?
      bool hasdKSq_;

      /// Will this WaveList be used for real-valued fields?
      bool isRealField_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      // Private member functions

      /// Access associated UnitCell<D> by reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Access associated Mesh<D> by reference.
      Mesh<D> const & mesh() const
      {  return *meshPtr_; }

   };

   // Get the array of minimum images on the device by reference.
   template <int D> inline
   DeviceArray<int> const & WaveList<D>::minImages_d() const
   {
      UTIL_CHECK(hasMinImages_);
      return minImages_;
   }

   // Get the kSq array on the device by reference.
   template <int D> inline
   RField<D> const & WaveList<D>::kSq() const
   {
      UTIL_CHECK(hasKSq_);
      return kSq_;
   }

   // Get a slice of the dKSq array on the device by reference.
   template <int D> inline 
   RField<D> const & WaveList<D>::dKSq(int i) const
   {
      UTIL_CHECK(hasdKSq_);
      return dKSqSlices_[i];
   }

   // Get the implicitInverse array by reference.
   template <int D> inline
   DeviceArray<bool> const & WaveList<D>::implicitInverse() const
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(isRealField_);
      return implicitInverse_;
   }

   #ifndef PRDC_CUDA_WAVE_LIST_TPP
   // Explicit instantiation declarations
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;
   #endif

} // Cuda
} // Prdc
} // Pscf
#endif

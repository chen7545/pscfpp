#ifndef RPC_CL_FIELD_IO_TPP
#define RPC_CL_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <prdc/cl/fieldIoUtil.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Log.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>

namespace Pscf {
namespace Prdc {
namespace Cl {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class CFT, class FFT>
   FieldIo<D,CFT,FFT>::FieldIo()
    : meshPtr_(nullptr),
      fftPtr_(nullptr),
      fileMasterPtr_(nullptr),
      nMonomer_(0),
      isAllocated_(false),
   {}

   /*
   * Destructor.
   */
   template <int D, class CFT, class FFT>
   FieldIo<D,CFT,FFT>::~FieldIo()
   {}

   // Initialization functions

   /*
   * Create associations with other members of a parent Domain<D> object.
   */
   template <int D, class CFT, class FFT>
   void
   FieldIo<D,CFT,FFT>::associate(
                    Mesh<D> const & mesh,
                    FFT const & fft,
                    typename UnitCell<D>::LatticeSystem const & lattice)
   {
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      latticePtr_ = &lattice;
   }

   /*
   * Create an association with a FileMaster.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::setFileMaster(
                              FileMaster const & fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Set nMonomer, the number of monomer types.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::setNMonomer(int nMonomer)
   {
      // Preconditions - require that function is only called once
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);

      nMonomer_ = nMonomer; 
   }

   /*
   * Open-close a file, and write an array of fields in r-grid format
   */
   template <int D, class CFT, class FFT>
   bool FieldIo<D,CFT,FFT>::readFields(
                              std::string filename,
                              DArray< CFT >& fields,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFields(file, fields, unitCell);
      file.close();
      return isSymmetric;
   }

   /*
   * Open-close a file, and read a single field in r-grid format
   */
   template <int D, class CFT, class FFT>
   bool FieldIo<D,CFT,FFT>::readField(
                              std::string filename,
                              CFT & field,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readField(file, field, unitCell);
      file.close();
      return isSymmetric;
   }

   /*
   * Open-close a file, and write an array of fields in r-grid format
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::writeFields(
                              std::string filename,
                              DArray< CFT > const & fields,
                              UnitCell<D> const & unitCell) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      bool writeHeader = true;
      bool writeMeshSize = true;
      writeFields(file, fields, unitCell,
                       writeHeader, writeMeshSize);
      file.close();
   }

   /*
   * Open-close a file, and write a single field in r-grid format
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::writeField(
                              std::string filename,
                              CFT const & field,
                              UnitCell<D> const & unitCell) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeField(file, field, unitCell, isSymmetric);
      file.close();
   }

   #if 0
   // FFT conversion functions (KGrid <-> RGrid) [Same in Rpc and Rpg]

   /*
   * Apply inverse FFT to an array of k-grid fields, converting to r-grid.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::convertKGridToRGrid(
                              DArray< CFT > const & in,
                              DArray< CFT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransform(in[i], out[i]);
      }
   }

   /*
   * Apply inverse FFT to a single k-grid field, converting to r-grid.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::convertKGridToRGrid(
                              CFT const & in, 
                              CFT& out) const
   {
      fft().inverseTransform(in, out);
   }

   /*
   * Apply forward FFT to an array of r-grid fields, converting to k-grid.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::convertRGridToKGrid(
                              DArray< CFT > const & in,
                              DArray< CFT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   /*
   * Apply forward FFT to a single r-grid field, converting to k-grid.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::convertRGridToKGrid(
                              CFT const & in,
                              CFT& out) const
   {  fft().forwardTransform(in, out); }

   /*
   * Convert a field file from k-grid to r-grid format.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::convertKGridToRGrid(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocate();
      UnitCell<D> tmpUnitCell;
      readFields(inFileName, tmpFields_, tmpUnitCell);
      for (int i = 0; i < nMonomer_; ++i) {
         fft().inverseTransform(tmpFields_[i],
                                tmpFields_[i]);
      }
      writeFields(outFileName, tmpFields_, tmpUnitCell);
   }

   /*
   * Convert a field file from r-grid to k-grid format.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::convertRGridToKGrid(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocate();
      UnitCell<D> tmpUnitCell;
      readFields(inFileName, tmpFields_, tmpUnitCell);
      convertRGridToKGrid(tmpFields_, tmpFieldsKGrid_);
      writeFieldsKGrid(outFileName, tmpFieldsKGrid_, tmpUnitCell);
   }

   // Field Inspection

   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::compareFields(
                                std::string const & filename1,
                                std::string const & filename2) const
   {
      DArray< CFT > fields1, fields2;
      UnitCell<D> tmpUnitCell;
      // Unallocated arrays will be allocated in readFields
      readFields(filename1, fields1, tmpUnitCell);
      readFields(filename2, fields2, tmpUnitCell);
      compareFields(fields1, fields2);
   }
   #endif

   // File Header IO Utilities

   /*
   * Read common part of field header.
   *
   * Extracts number of monomers (i.e., number of fields) and unitCell
   * data (lattice type and parameters) from the file, and returns these
   * via non-const reference parameters nMonomer and unitCell.
   *
   * Also validates lattice type and groupName (if any), and constructs
   * a symmetry-adapted basis if there is a group but no basis.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::readFieldHeader(
                              std::istream& in,
                              int& nMonomer,
                              UnitCell<D>& unitCell) const
   {
      // Preconditions
      UTIL_CHECK(latticePtr_);
      if (unitCell.lattice() == UnitCell<D>::Null) {
         UTIL_CHECK(unitCell.nParameter() == 0);
      } else {
         UTIL_CHECK(unitCell.nParameter() > 0);
         UTIL_CHECK(unitCell.lattice() == lattice());
      }

      // Read field header to set unitCell, groupNameIn, nMonomer
      int ver1, ver2;
      std::string groupNameIn;

      Pscf::Prdc::readFieldHeader(in, ver1, ver2, unitCell,
                                  groupNameIn, nMonomer);
      // Note: Function definition in prdc/crystal/fieldHeader.tpp

      // Checks of data from header
      UTIL_CHECK(ver1 == 1);
      //UTIL_CHECK(ver2 == 0);
      UTIL_CHECK(unitCell.isInitialized());
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.nParameter() > 0);

      // Validate or initialize lattice type
      if (lattice() != unitCell.lattice()) {
         Log::file() << std::endl
               << "Error - Mismatched lattice types "
               << "in FieldIo function readFieldHeader:\n"
               << "  FieldIo::lattice  :" << lattice() << "\n"
               << "  Unit cell lattice :" << unitCell.lattice()
               << "\n";
         UTIL_THROW("Mismatched lattice types");
      }

   }

   /*
   * Write a field file header.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::writeFieldHeader(
                              std::ostream &out,
                              int nMonomer,
                              UnitCell<D> const & unitCell) const
   {
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      if (isSymmetric) {
         UTIL_CHECK(hasGroup());
         gName = groupName();
      }
      Pscf::Prdc::writeFieldHeader(out, v1, v2, unitCell,
                                   gName, nMonomer);
      // Note: This function is defined in prdc/crystal/fieldHeader.tpp
   }

   // Protected functions to check and allocate private workspace arrays

   /*
   * If necessary, allocate r-grid workspace.
   */
   template <int D, class CFT, class FFT>
   void FieldIo<D,CFT,FFT>::checkAllocate() const
   {
      if (isAllocated_) return;

      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(mesh().size() > 0);
      IntVec<D> const & meshDimensions = mesh().dimensions();
      tmpFields_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         tmpFields_[i].allocate(meshDimensions);
      }
      isAllocated_ = true;
   }

} // namespace Cl
} // namespace Prdc
} // namespace Pscf
#endif

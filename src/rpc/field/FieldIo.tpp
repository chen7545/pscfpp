#ifndef RPC_FIELD_IO_TPP
#define RPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/crystal/fieldHeader.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/misc/FileMaster.h>
#include <util/misc/Log.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : meshPtr_(0),
      fftPtr_(0),
      hasGroupPtr_(0),
      groupNamePtr_(0),
      groupPtr_(0),
      basisPtr_(0),
      fileMasterPtr_()
   {}

   /*
   * Destructor.
   */
   template <int D>
   FieldIo<D>::~FieldIo()
   {}

   /*
   * Create associations with other members of parent Domain.
   */
   template <int D>
   void
   FieldIo<D>::associate(
                    Mesh<D> const & mesh,
                    FFT<D> const & fft,
                    typename UnitCell<D>::LatticeSystem const & lattice,
                    bool const & hasGroup,
                    std::string const & groupName,
                    SpaceGroup<D> const & group,
                    Basis<D> & basis)
   {
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      latticePtr_ = &lattice;
      hasGroupPtr_ = &hasGroup;
      groupNamePtr_ = &groupName;
      groupPtr_ = &group;
      basisPtr_ = &basis;
   }

   /*
   * Create an association with a FileMaster.
   */
   template <int D>
   void FieldIo<D>::setFileMaster(FileMaster const & fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Read a set of fields in basis format.
   */
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::istream& in,
                                    DArray< DArray<double> >& fields,
                                    UnitCell<D>& unitCell) const
   {
      // Precondition
      UTIL_CHECK(hasGroup());

      // Read header of field file
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(isSymmetric);
      UTIL_CHECK(basis().isInitialized());

      // Read the number of stars into nStarIn
      std::string label;
      in >> label;
      if (label != "N_star" && label != "N_basis") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected N_basis or N_star, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      int nStarIn;
      in >> nStarIn;
      UTIL_CHECK(nStarIn > 0);

      // If "fields" parameter is allocated, check if dimensions
      // match those of the system's mesh.  Otherwise, allocate.

      int fieldCapacity;
      if (fields.isAllocated()) {
         // If outer DArray is allocated, require that it matches the
         // number of inputted fields and that internal DArrays are also
         // allocated all with the same dimensions.
         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0);
         UTIL_CHECK(nMonomerFields == nMonomer);

         // Check that all internal DArrays have same dimension.
         fieldCapacity = fields[0].capacity();
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK( fields[i].capacity() == fieldCapacity );
         }
      } else {
         // Else, allocate fields to the number of inputted fields
         // and the internal dimensions to the number of stars.
         fields.allocate(nMonomer);
         fieldCapacity = nStarIn;

         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(fieldCapacity);
         }
      }

      // Initialize all field array elements to zero
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < fieldCapacity; ++j) {
            fields[i][j] = 0.0;
         }
      }

      // Reset nStarIn = min(nStarIn, fieldCapacity)
      if (fieldCapacity < nStarIn) {
         nStarIn = fieldCapacity;
      }

      // Allocate temp arrays used to read in components
      DArray<double> temp, temp2;
      temp.allocate(nMonomer);
      temp2.allocate(nMonomer);

      typename Basis<D>::Star const * starPtr;
      typename Basis<D>::Star const * starPtr2;
      IntVec<D> waveIn, waveIn2;
      int sizeIn, sizeIn2;
      int starId, starId2;
      int basisId, basisId2;
      int waveId, waveId2;

      std::complex<double> coeff, phasor;
      IntVec<D> waveBz, waveDft;
      int nReversedPair = 0;
      bool waveExists, sizeMatches;

      // Loop over stars in input file to read field components
      int i = 0;
      while (i < nStarIn) {

         // Read next line of data
         for (int j = 0; j < nMonomer; ++j) {
            in >> temp[j];        // field components
         }
         in >> waveIn;            // wave of star
         in >> sizeIn;            // # of waves in star
         ++i;

         sizeMatches = false;
         waveExists = false;

         // Check if waveIn is in first Brillouin zone (FBZ) for the mesh.
         waveBz = shiftToMinimum(waveIn, mesh().dimensions(), unitCell);
         waveExists = (waveIn == waveBz);

         if (waveExists) {

            // Find the star containing waveIn
            waveDft = waveIn;
            mesh().shift(waveDft);
            waveId = basis().waveId(waveDft);
            starId = basis().wave(waveId).starId;
            starPtr = &basis().star(starId);
            UTIL_CHECK(!(starPtr->cancel));
            basisId = starPtr->basisId;

            if (starPtr->size == sizeIn) {
               sizeMatches = true;
            } else {
               Log::file()
                  <<  "Warning: Inconsistent star size (line ignored)\n"
                  <<  "wave from file = " << waveIn << "\n"
                  <<  "size from file = " << sizeIn << "\n"
                  <<  "size of star   = " << starPtr->size
                  << "\n";
               sizeMatches = false;
            }

         }

         if (waveExists && sizeMatches) { // Attempt to process wave

            if (starPtr->invertFlag == 0) {

               if (starPtr->waveBz == waveIn) {

                  // Copy components of closed star to fields array
                  for (int j = 0; j < nMonomer; ++j) {
                      fields[j][basisId] = temp[j];
                  }

               } else {

                  Log::file()
                     <<  "Inconsistent wave of closed star on input\n"
                     <<  "wave from file = " << waveIn  << "\n"
                     <<  "starId of wave = " << starId  << "\n"
                     <<  "waveBz of star = " << starPtr->waveBz
                     << "\n";

               }

            } else {

               // Read the next line
               for (int j = 0; j < nMonomer; ++j) {
                  in >> temp2[j];          // components of field
               }
               in >> waveIn2;              // wave of star
               in >> sizeIn2;              // # of wavevectors in star
               ++i;

               // Check that waveIn2 is also in the 1st BZ
               waveBz =
                  shiftToMinimum(waveIn2, mesh().dimensions(), unitCell);
               UTIL_CHECK(waveIn2 == waveBz);

               // Identify the star containing waveIn2
               waveDft = waveIn2;
               mesh().shift(waveDft);
               waveId2 = basis().waveId(waveDft);
               starId2 = basis().wave(waveId2).starId;
               starPtr2 = &basis().star(starId2);
               UTIL_CHECK(!(starPtr2->cancel));
               basisId2 = starPtr2->basisId;
               UTIL_CHECK(starPtr2->size == sizeIn2);
               UTIL_CHECK(sizeIn == sizeIn2);

               if (starPtr->invertFlag == 1) {

                  // This is a pair of open stars written in the same
                  // order as in this basis. Check preconditions:
                  UTIL_CHECK(starPtr2->invertFlag == -1);
                  UTIL_CHECK(starId2 = starId + 1);
                  UTIL_CHECK(basisId2 = basisId + 1);
                  UTIL_CHECK(starPtr->waveBz == waveIn);
                  UTIL_CHECK(starPtr2->waveBz == waveIn2);

                  // Copy components for both stars into fields array
                  for (int j = 0; j < nMonomer; ++j) {
                      fields[j][basisId] = temp[j];
                      fields[j][basisId2] = temp2[j];
                  }

               } else
               if (starPtr->invertFlag == -1) {

                  // This is a pair of open stars written in opposite
                  // order from in this basis. Check preconditions:
                  UTIL_CHECK(starPtr2->invertFlag == 1);
                  UTIL_CHECK(starId == starId2 + 1);
                  UTIL_CHECK(basisId == basisId2 + 1);
                  UTIL_CHECK(waveId == starPtr->beginId);

                  // Check that waveIn2 is negation of waveIn
                  IntVec<D> nVec;
                  nVec.negate(waveIn);
                  nVec =
                       shiftToMinimum(nVec, mesh().dimensions(), unitCell);
                  UTIL_CHECK(waveIn2 == nVec);

                  /*
                  * Consider two related stars, C and D, that are listed
                  * in the order (C,D) in the basis used in this code (the
                  * reading program), but that were listed in the opposite
                  * order (D,C) in the program that wrote the file (the
                  * writing program). In the basis of the reading program,
                  * star C has star index starId2, while star D has index
                  * starId = starid2 + 1.
                  *
                  * Let f(r) and f^{*}(r) denote the basis functions used
                  * by the reading program for stars C and D, respectively.
                  * Let u(r) and u^{*}(r) denote the corresponding basis
                  * functions used by the writing program for stars C
                  * and D.  Let exp(i phi) denote the unit magnitude
                  * coefficient (i.e., phasor) within f(r) of the wave
                  * with wave index waveId2, which was the characteristic
                  * wave for star C in the writing program. The
                  * coefficient of this wave within the basis function
                  * u(r) used by the writing program must instead be real
                  * and positive. This implies that
                  * u(r) = exp(-i phi) f(r).
                  *
                  * Focus on the contribution to the field for a specific
                  * monomer type j.  Let a and b denote the desired
                  * coefficients of stars C and D in the reading program,
                  * for which the total contribution of both stars to the
                  * field is:
                  *
                  *  (a - ib) f(r) + (a + ib) f^{*}(r)
                  *
                  * Let A = temp[j] and B = temp2[j] denote the
                  * coefficients read from file in order (A,B).  Noting
                  * that the stars were listed in the order (D,C) in the
                  * basis used by the writing program, the contribution
                  * of both stars must be (A-iB)u^{*}(r)+(A+iB)u(r), or:
                  *
                  *  (A+iB) exp(-i phi)f(r) + (A-iB) exp(i phi) f^{*}(r)
                  *
                  * Comparing coefficients of f^{*}(r), we find that
                  *
                  *       (a + ib) = (A - iB) exp(i phi)
                  *
                  * This equality is implemented below, where the
                  * variable "phasor" is set equal to exp(i phi).
                  */
                  phasor = basis().wave(waveId2).coeff;
                  phasor = phasor/std::abs(phasor);
                  for (int j = 0; j < nMonomer; ++j) {
                      coeff = std::complex<double>(temp[j],-temp2[j]);
                      coeff *= phasor;
                      fields[j][basisId2] = real(coeff);
                      fields[j][basisId ] = imag(coeff);
                  }

                  // Increment count of number of reversed open pairs
                  ++nReversedPair;

               } else {
                  UTIL_THROW("Invalid starInvert value");
               }

            }   // if (wavePtr->invertFlag == 0) ... else ...

         }   // if (waveExists && sizeMatches)

      }   // end while (i < nStarIn)

      if (nReversedPair > 0) {
         Log::file() << "\n";
         Log::file() << nReversedPair << " reversed pairs of open stars"
                   << " detected in FieldIo::readFieldsBasis\n";
      }

   }

   /*
   * Open/close a file and read a set of fields in basis format.
   */
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::string filename,
                                    DArray<DArray<double> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {

      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsBasis(file, fields, unitCell);
      file.close();
   }

   /*
   * Read a single fields in basis format into input stream
   */
   template <int D>
   void FieldIo<D>::readFieldBasis(std::istream& in,
                                   DArray<double>& field,
                                   UnitCell<D>& unitCell) const
   {
      // Local array, of data type required by readFieldsBasis
      DArray<DArray<double> > fields;

      // If field is allocated, allocate local array fields
      // Otherwise, pass unallocated fields array to readFieldsBasis
      if (field.isAllocated()) {
         fields.allocate(1);
         fields[0].allocate(field.capacity());
      }

      // Read file containing a single field, allocate if needed.
      readFieldsBasis(in, fields, unitCell);
      // Check that it only read 1 field
      UTIL_CHECK(fields.capacity() == 1);

      // Copy data from local array to function parameter
      field = fields[0];
   }

   /*
   * Open-close a file and read single fields in basis format.
   */
   template <int D>
   void FieldIo<D>::readFieldBasis(std::string filename,
                                   DArray<double>& field,
                                   UnitCell<D>& unitCell)
   const
   {
      // Local array, of data type required by readFieldsBasis
      DArray<DArray<double> > fields;

      // If field is allocated, allocate local array fields
      // Otherwise, pass unallocated fields array to readFieldsBasis
      if (field.isAllocated()) {
         fields.allocate(1);
         fields[0].allocate(field.capacity());
      }

      // Read file containing a single field, allocate if needed.
      readFieldsBasis(filename, fields, unitCell);
      // Check that it only read 1 field
      UTIL_CHECK(fields.capacity() == 1);

      // Copy data from local array to function parameter
      field = fields[0];
   }

   /*
   * Write an array of fields in basis format to an output stream.
   */
   template <int D>
   void
   FieldIo<D>::writeFieldsBasis(std::ostream &out,
                                DArray<DArray<double> > const & fields,
                                UnitCell<D> const & unitCell) const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(basis().isInitialized());

      // Write header (common portion)
      bool isSymmetric = true;
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);

      // Write N_Basis
      int nStar = basis().nStar();
      int nBasis = basis().nBasis();
      out << "N_basis      " << std::endl
          << "             " << nBasis << std::endl;

      // Write fields
      int ib = 0;
      for (int i = 0; i < nStar; ++i) {
         if (!basis().star(i).cancel) {
            for (int j = 0; j < nMonomer; ++j) {
               out << Dbl(fields[j][ib], 20, 10);
            }
            out << "   ";
            for (int j = 0; j < D; ++j) {
               out << Int(basis().star(i).waveBz[j], 5);
            }
            out << Int(basis().star(i).size, 5) << std::endl;
            ++ib;
         }
      }

   }

   /*
   * Open-close a file, and write an array of fields in basis format.
   */
   template <int D>
   void
   FieldIo<D>::writeFieldsBasis(std::string filename,
                                DArray<DArray<double> > const & fields,
                                UnitCell<D> const & unitCell)
   const
   {
       std::ofstream file;
       fileMaster().openOutputFile(filename, file);
       writeFieldsBasis(file, fields, unitCell);
       file.close();
   }

   /*
   * Write a single field in basis format to an output stream.
   */
   template <int D>
   void FieldIo<D>::writeFieldBasis(std::ostream& out,
                                    DArray<double> const & field,
                                    UnitCell<D> const & unitCell)
   const
   {
      // Create local array of type required by writeFieldsBasis
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());

      // Copy data from input parameter to local array
      fields[0] = field;

      writeFieldsBasis(out, fields, unitCell);
   }

   template <int D>
   void
   FieldIo<D>::writeFieldBasis(std::string filename,
                               DArray<double> const & field,
                               UnitCell<D> const & unitCell) const
   {
      // Create local array of type required by writeFieldsBasis
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());

      // Copy data from input parameter to local array
      fields[0] = field;

      writeFieldsBasis(filename, fields, unitCell);
   }

   // R-Grid Field Format IO

   /*
   * Read an array of fields in r-grid format from an input stream.
   */
   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::istream &in,
                                    DArray<RField<D> >& fields,
                                    UnitCell<D>& unitCell) const
   {
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);

      // If "fields" parameter is allocated, check if dimensions match
      // those of the system's mesh.  Otherwise, allocate.

      if (fields.isAllocated()) {
         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0);
         UTIL_CHECK(nMonomerFields == nMonomer);

         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].meshDimensions()==mesh().dimensions());
         }
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(mesh().dimensions());
         }
      }

      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> nGrid;
      in >> nGrid;
      if (nGrid != mesh().dimensions()) {
         Log::file()
             << "In Rpc::FieldIo::readFieldsRGrid:\n"
             << "mesh().dimensions() = " << mesh().dimensions() << "\n"
             << "nGrid               = " << mesh().dimensions() << "\n";
         UTIL_THROW("Inconsistent grid in mesh and field header");
      }

      // Setup temporary workspace array.
      DArray<RField<D> > temp;
      temp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         temp[i].allocate(mesh().dimensions());
      }

      // Read Fields;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nMonomer; ++i) {
            in  >> std::setprecision(15) >> temp[i][itr.rank()];
         }
      }

      int p = 0;
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;

      if (D==3) {
         while (n1 < mesh().dimension(0)) {
            q = p;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               r = q;
               n3 = 0;
               while (n3 < mesh().dimension(2)) {
                  for (int i = 0; i < nMonomer; ++i) {
                     fields[i][s] = temp[i][r];
                  }
                  r = r + (mesh().dimension(0) * mesh().dimension(1));
                  ++s;
                  ++n3;
               }
               q = q + mesh().dimension(0);
               ++n2;
            }
            ++n1;
            ++p;
         }
      }

      else if (D==2) {
         while (n1 < mesh().dimension(0)) {
            r =q;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               for (int i = 0; i < nMonomer; ++i) {
                  fields[i][s] = temp[i][r];
               }
               r = r + (mesh().dimension(0));
               ++s;
               ++n2;
            }
            ++q;
            ++n1;
         }
      }

      else if (D==1) {

         while (n1 < mesh().dimension(0)) {
            for (int i = 0; i < nMonomer; ++i) {
               fields[i][s] = temp[i][r];
            }
            ++r;
            ++s;
            ++n1;
         }
      }

      else{
         Log::file() << "Invalid Dimensions";
      }

   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::string filename,
                              DArray< RField<D> >& fields,
                              UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsRGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::istream &in,
                                    RField<D> & field,
                                    UnitCell<D>& unitCell)
   const
   {
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);

      // Only reading in a file with a single field.
      UTIL_CHECK(nMonomer == 1);

      // If "field" parameter is allocated, check if dimensions match
      // those of the system's mesh.  Otherwise, allocate.

      if (field.isAllocated()) {
         UTIL_CHECK(field.meshDimensions() == mesh().dimensions());
      } else {
         field.allocate(mesh().dimensions());
      }

      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      // Setup temporary workspace.
      RField<D> temp;
      temp.allocate(mesh().dimensions());

      // Read Field;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         in  >> std::setprecision(15) >> temp[itr.rank()];
      }

      int p = 0;
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;

      if (D==3) {
         while (n1 < mesh().dimension(0)) {
            q = p;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               r = q;
               n3 = 0;
               while (n3 < mesh().dimension(2)) {
                  field[s] = temp[r];
                  r = r + (mesh().dimension(0) * mesh().dimension(1));
                  ++s;
                  ++n3;
               }
               q = q + mesh().dimension(0);
               ++n2;
            }
            ++n1;
            ++p;
         }
      }

      else if (D==2) {
         while (n1 < mesh().dimension(0)) {
            r =q;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               field[s] = temp[r];
               r = r + (mesh().dimension(0));
               ++s;
               ++n2;
            }
            ++q;
            ++n1;
         }
      }

      else if (D==1) {

         while (n1 < mesh().dimension(0)) {
               field[s] = temp[r];
            ++r;
            ++s;
            ++n1;
         }
      }

      else{
         Log::file() << "Invalid Dimensions";
      }
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::string filename,
                                    RField<D> & field,
                                    UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldRGrid(file, field, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsRGridData(std::istream& in,
                                        DArray< RField<D> >& fields,
                                        int nMonomer)
   const
   {
      // If "fields" parameter is allocated, check if mesh size match
      // those of the system's mesh.  Otherwise, allocate.

      if (fields.isAllocated()) {
         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0);
         UTIL_CHECK(nMonomerFields == nMonomer);

         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].capacity() == mesh().size());
         }
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(mesh().dimensions());
         }
      }

      // Setup temporary workspace array.
      DArray<RField<D> > temp;
      temp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         temp[i].allocate(mesh().dimensions());
      }

      // Read Fields;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nMonomer; ++i) {
            in  >> std::setprecision(15) >> temp[i][itr.rank()];
         }
      }

      int p = 0;
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;

      if (D==3) {
         while (n1 < mesh().dimension(0)) {
            q = p;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               r = q;
               n3 = 0;
               while (n3 < mesh().dimension(2)) {
                  for (int i = 0; i < nMonomer; ++i) {
                     fields[i][s] = temp[i][r];
                  }
                  r = r + (mesh().dimension(0) * mesh().dimension(1));
                  ++s;
                  ++n3;
               }
               q = q + mesh().dimension(0);
               ++n2;
            }
            ++n1;
            ++p;
         }
      }

      else if (D==2) {
         while (n1 < mesh().dimension(0)) {
            r =q;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               for (int i = 0; i < nMonomer; ++i) {
                  fields[i][s] = temp[i][r];
               }
               r = r + (mesh().dimension(0));
               ++s;
               ++n2;
            }
            ++q;
            ++n1;
         }
      }

      else if (D==1) {

         while (n1 < mesh().dimension(0)) {
            for (int i = 0; i < nMonomer; ++i) {
               fields[i][s] = temp[i][r];
            }
            ++r;
            ++s;
            ++n1;
         }
      }

      else{
         Log::file() << "Invalid Dimensions";
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::ostream &out,
                                     DArray<RField<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     bool writeHeader,
                                     bool isSymmetric,
                                     bool writeMeshSize) const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      if (writeHeader){
         writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      }
      
      IntVec<D> meshDimensions = fields[0].meshDimensions();
      if (writeMeshSize){
         out << "mesh " <<  std::endl
             << "           " << meshDimensions << std::endl;
      }
      
      DArray<RField<D> > temp;
      temp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {

         temp[i].allocate(meshDimensions);
      } 

      int p = 0;
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 =0;
      int n2 =0;
      int n3 =0;

      if (D==3) {
         while (n3 < meshDimensions[2]) {
            q = p; 
            n2 = 0; 
            while (n2 < meshDimensions[1]) {
               r =q;
               n1 = 0; 
               while (n1 < meshDimensions[0]) {
                  for (int i = 0; i < nMonomer; ++i) {
                     temp[i][s] = fields[i][r];
                  }    
                  r = r + (meshDimensions[1] * meshDimensions[2]);
                  ++s; 
                  ++n1;     
               }    
               q = q + meshDimensions[2];

               ++n2;
            }
            ++n3;
            ++p;
         }
      }
      else if (D==2) {
         while (n2 < meshDimensions[1]) {
            r =q;
            n1 = 0;
            while (n1 < meshDimensions[0]) {
               for (int i = 0; i < nMonomer; ++i) {
                  temp[i][s] = fields[i][r];
               }
               r = r + (meshDimensions[1]);
               ++s;
               ++n1;
            }
            ++q;
            ++n2;
         }
      }
      else if (D==1) {
         while (n1 < meshDimensions[0]) {
            for (int i = 0; i < nMonomer; ++i) {
               temp[i][s] = fields[i][r];
            }
            ++r;
            ++s;
            ++n1;
         }
      } else {
         Log::file() << "Invalid Dimensions";
      }

      // Write fields
      MeshIterator<D> itr(meshDimensions);
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int j = 0; j < nMonomer; ++j) {
            out << "  " << Dbl(temp[j][itr.rank()], 18, 15);
         }
         out << std::endl;
      }

   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::string filename,
                                     DArray< RField<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      bool writeHeader = true;
      bool writeMeshSize = true;
      writeFieldsRGrid(file, fields, unitCell, writeHeader, isSymmetric, writeMeshSize);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::ostream &out,
                                    RField<D> const & field,
                                    UnitCell<D> const & unitCell,
                                    bool writeHeader,
                                    bool isSymmetric)
   const
   {
      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell, isSymmetric);
         out << "mesh " <<  std::endl
             << "           " << mesh().dimensions() << std::endl;
      }

      RField<D> temp;
      temp.allocate(mesh().dimensions());

      int p = 0;
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 =0;
      int n2 =0;
      int n3 =0;

      if (D==3) {
         while (n3 < mesh().dimension(2)) {
            q = p;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               r =q;
               n1 = 0;
               while (n1 < mesh().dimension(0)) {
                  temp[s] = field[r];
                  r = r + (mesh().dimension(1) * mesh().dimension(2));
                  ++s;
                  ++n1;
               }
               q = q + mesh().dimension(2);
               ++n2;
            }
            ++n3;
            ++p;
         }
      }
      else if (D==2) {
         while (n2 < mesh().dimension(1)) {
            r =q;
            n1 = 0;
            while (n1 < mesh().dimension(0)) {
               temp[s] = field[r];
               r = r + (mesh().dimension(1));
               ++s;
               ++n1;
            }
            ++q;
            ++n2;
         }
      }
      else if (D==1) {
         while (n1 < mesh().dimension(0)) {
            temp[s] = field[r];
            ++r;
            ++s;
            ++n1;
         }
      } else {
         Log::file() << "Invalid Dimensions";
      }

      // Write field
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         out << "  " << Dbl(temp[itr.rank()], 18, 15);
         out << std::endl;
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::string filename,
                                    RField<D> const & field,
                                    UnitCell<D> const & unitCell,
                                    bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldRGrid(file, field, unitCell, isSymmetric);
      file.close();
   }

   // K-Grid Field Format

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::istream &in,
                                    DArray<RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);

      // If "fields" parameter is allocated, check if dimensions match
      // those of the system's mesh.  Otherwise, allocate.

      if (fields.isAllocated()) {

         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0)
         UTIL_CHECK(nMonomerFields == nMonomer)

         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
         }

      } else {

         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(mesh().dimensions());
         }

      }

      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      // Read fields;
      int i, idum;
      MeshIterator<D> itr(fields[0].dftDimensions());
      i = 0;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         in >> idum;
         UTIL_CHECK(i == idum);
         UTIL_CHECK(i == itr.rank());
         for (int i = 0; i < nMonomer; ++i) {
            for (int j = 0; j < 2; ++j) {
               in >> fields[i][itr.rank()][j];
            }
         }
         ++i;
      }
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::string filename,
                                    DArray< RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields, unitCell);
      file.close();
   }
   
   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::ostream &out,
                                     DArray<RFieldDft<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     bool isSymmetric) const
   {
      // Inspect fields array
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
      }

      // Write header
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      out << "mesh " << std::endl
          << "               " << mesh().dimensions() << std::endl;

      // Write fields
      MeshIterator<D> itr(fields[0].dftDimensions());
      int i = 0;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         UTIL_CHECK(i == itr.rank());
         out << Int(itr.rank(), 5);
         for (int j = 0; j < nMonomer; ++j) {
               out << "  "
                   << Dbl(fields[j][itr.rank()][0], 20, 12)
                   << Dbl(fields[j][itr.rank()][1], 20, 12);
         }
         out << std::endl;
         ++i;
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::string filename,
                                    DArray< RFieldDft<D> > const & fields,
                                    UnitCell<D> const & unitCell,
                                    bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields, unitCell, isSymmetric);
      file.close();
   }

   /*
   * Read common part of field header and extract
   * the number of monomers (number of fields) in the file.
   */
   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in,
                                    int& nMonomer,
                                    UnitCell<D>& unitCell,
                                    bool & isSymmetric)
   const
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
               << "in function FieldIo<D>::readFieldHeader:\n"
               << "  FieldIo::lattice  :" << lattice() << "\n"
               << "  Unit cell lattice :" << unitCell.lattice()
               << "\n";
         UTIL_THROW("Mismatched lattice types");
      }

      // Check for presence of group name
      isSymmetric = false;
      if (groupNameIn != "") {
         isSymmetric = true;
      }

      // Process group and basis (if any)
      if (hasGroup()) {

         // Check consistency of groupName values
         if (isSymmetric) {
            UTIL_CHECK(groupNamePtr_);
            if (groupNameIn != groupName()) {
               Log::file() << std::endl
                  << "Error - Mismatched group names in "
                  << "function FieldIo<D>::readFieldHeader:\n"
                  << "  FieldIo::groupName :" << groupName() << "\n"
                  << "  Field file header  :" << groupNameIn << "\n";
               UTIL_THROW("Mismatched group names");
            }
         }

         // If there is a group but no basis, construct a basis
         UTIL_CHECK(basisPtr_);
         if (!basis().isInitialized()) {
            basisPtr_->makeBasis(mesh(), unitCell, group());
         }
         UTIL_CHECK(basis().isInitialized());

      }

   }

   template <int D>
   void FieldIo<D>::writeFieldHeader(std::ostream &out,
                                     int nMonomer,
                                     UnitCell<D> const & unitCell,
                                     bool isSymmetric) const
   {
      int v1 = 1;
      int v2 = 0;
      std::string gname = "";
      if (isSymmetric) {
         UTIL_CHECK(hasGroup());
         gname = groupName();
      }
      Pscf::Prdc::writeFieldHeader(out, v1, v2, unitCell,
                                   gname, nMonomer);
      // Note: This function is defined in prdc/crystal/fieldHeader.tpp
   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray<double> const & in,
                                        RFieldDft<D>& out) const
   {
      UTIL_CHECK(basis().isInitialized());

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(out.dftDimensions());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient for star
      std::complex<double> coeff;             // coefficient for wave
      IntVec<D> indices;                      // dft grid indices of wave
      int rank;                               // dft grid rank of wave
      int is;                                 // star index
      int ib;                                 // basis index
      int iw;                                 // wave index

      // Initialize all dft coponents to zero
      for (rank = 0; rank < dftMesh.size(); ++rank) {
         out[rank][0] = 0.0;
         out[rank][1] = 0.0;
      }

      // Loop over stars, skipping cancelled stars
      is = 0;
      while (is < basis().nStar()) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         // Set basisId for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Make complex coefficient for star basis function
            component = std::complex<double>(in[ib], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!wavePtr->implicit) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Loop over waves in first star
            component = std::complex<double>(in[ib], -in[ib+1]);
            component /= sqrt(2.0);
            starPtr = &(basis().star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis().star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = std::complex<double>(in[ib], +in[ib+1]);
            component /= sqrt(2.0);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }

            // Increment is by 2 (two stars were processed)
            is += 2;

         } else {

            UTIL_THROW("Invalid invertFlag value");

         }

      }

   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(RFieldDft<D> const & in,
                                        DArray<double>& out,
                                        bool checkSymmetry,
                                        double epsilon) const
   {
      UTIL_CHECK(basis().isInitialized());

      if (checkSymmetry) {
         // Check if kgrid has symmetry
         bool symmetric = hasSymmetry(in, epsilon, true);
         if (!symmetric) {
            Log::file() << std::endl
               << "WARNING: non-negligible error in conversion to "
               << "symmetry-adapted basis format." << std::endl
               << "   See error values printed above for each "
               << "asymmetric field." << std::endl
               << "   The field that is output by the above operation "
               << "will be a" << std::endl
               << "   symmetrized version of the input field."
               << std::endl << std::endl;
         }
      }

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(in.dftDimensions());

      typename Basis<D>::Star const* starPtr;  // pointer to current star
      typename Basis<D>::Wave const* wavePtr;  // pointer to current wave
      std::complex<double> component;          // coefficient for star
      int rank;                                // dft grid rank of wave
      int is;                                  // star index
      int ib;                                  // basis index

      // Initialize all components to zero
      for (is = 0; is < basis().nBasis(); ++is) {
         out[is] = 0.0;
      }

      // Loop over stars
      is = 0;
      while (is < basis().nStar()) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         // Set basis id for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Choose a wave in the star that is not implicit
            int beginId = starPtr->beginId;
            int endId = starPtr->endId;
            int iw = 0;
            bool isImplicit = true;
            while (isImplicit) {
               wavePtr = &basis().wave(beginId + iw);
               if (!wavePtr->implicit) {
                  isImplicit = false;
               } else {
                   UTIL_CHECK(beginId + iw < endId - 1 - iw);
                   wavePtr = &basis().wave(endId - 1 - iw);
                   if (!wavePtr->implicit) {
                      isImplicit = false;
                   }
               }
               ++iw;
            }
            UTIL_CHECK(wavePtr->starId == is);

            // Compute component value
            rank = dftMesh.rank(wavePtr->indicesDft);
            component = std::complex<double>(in[rank][0], in[rank][1]);
            component /= wavePtr->coeff;
            out[ib] = component.real();
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Identify a characteristic wave that is not implicit:
            // Either the first wave of the 1st star or last wave of 2nd
            wavePtr = &basis().wave(starPtr->beginId);
            UTIL_CHECK(wavePtr->starId == is);
            if (wavePtr->implicit) {
               starPtr = &(basis().star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &basis().wave(starPtr->endId - 1);
               UTIL_CHECK(!(wavePtr->implicit));
               UTIL_CHECK(wavePtr->starId == is+1);
            }
            rank = dftMesh.rank(wavePtr->indicesDft);
            component = std::complex<double>(in[rank][0], in[rank][1]);
            UTIL_CHECK(std::abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);

            // Compute basis function coefficient values
            if (starPtr->invertFlag == 1) {
               out[ib] = component.real();
               out[ib+1] = -component.imag();
            } else {
               out[ib] = component.real();
               out[ib+1] = component.imag();
            }

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is
   }

   template <int D>
   void
   FieldIo<D>::convertBasisToKGrid(DArray< DArray <double> > const & in,
                                   DArray< RFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }
   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(DArray< RFieldDft<D> > const & in,
                                        DArray< DArray <double> > & out,
                                        bool checkSymmetry,
                                        double epsilon) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();

      bool symmetric(true);
      for (int i = 0; i < n; ++i) {
         if (checkSymmetry) {
            // Check if kgrid has symmetry
            bool tmp_sym = hasSymmetry(in[i], epsilon, true);
            if (!tmp_sym) symmetric = false;
         }
         convertKGridToBasis(in[i], out[i], false);
      }

      // Print warning if any input field is assymmetric
      if (!symmetric) {
         Log::file() << std::endl
            << "WARNING: non-negligible error in conversion to "
            << "symmetry-adapted basis format." << std::endl
            << "See error values printed above for each asymmetric field."
            << std::endl
            << "The field that is output by this operation will be "
            << "a symmetrized version of" << std::endl
            << "the input field." << std::endl << std::endl;
      }
   }

   template <int D>
   void
   FieldIo<D>::convertBasisToRGrid(DArray<double> const & in,
                                   RField<D>& out) const
   {
      checkWorkDft();
      convertBasisToKGrid(in, workDft_);
      fft().inverseTransformSafe(workDft_, out);
   }

   template <int D>
   void
   FieldIo<D>::convertBasisToRGrid(DArray< DArray <double> > const & in,
                                   DArray< RField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkWorkDft();

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], workDft_);
         fft().inverseTransformSafe(workDft_, out[i]);
      }
   }

   template <int D>
   void
   FieldIo<D>::convertRGridToBasis(RField<D> const & in,
                                   DArray<double> & out,
                                   bool checkSymmetry,
                                   double epsilon) const
   {
      checkWorkDft();
      fft().forwardTransform(in, workDft_);
      convertKGridToBasis(workDft_, out, checkSymmetry, epsilon);
   }

   template <int D>
   void
   FieldIo<D>::convertRGridToBasis(DArray< RField<D> > const & in,
                                   DArray< DArray <double> > & out,
                                   bool checkSymmetry,
                                   double epsilon) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkWorkDft();

      int n = in.capacity();

      bool symmetric(true);
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], workDft_);
         if (checkSymmetry) {
            // Check if kgrid has symmetry
            bool tmp_sym = hasSymmetry(workDft_, epsilon, true);
            if (!tmp_sym) symmetric = false;
         }
         convertKGridToBasis(workDft_, out[i], false);
      }

      // Print warning if any input fields is asymmetric
      if (!symmetric) {
         Log::file() << std::endl
             << "WARNING: non-negligible error in conversion to "
             << "symmetry-adapted basis format." << std::endl
             << "   See error values printed above for each "
             << "asymmetric field." << std::endl
             << "   The field that is output by the above operation "
             << "will be a" << std::endl
             << "   symmetrized version of the input field."
             << std::endl << std::endl;
      }
   }

   /*
   * Apply inverse FFT to an array of k-grid fields.
   */
   template <int D>
   void
   FieldIo<D>::convertKGridToRGrid(DArray< RFieldDft<D> > & in,
                                   DArray< RField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransformSafe(in[i], out[i]);
      }
   }

   /*
   * Apply inverse FFT to a single k-grid field.
   */
   template <int D>
   void
   FieldIo<D>::convertKGridToRGrid(RFieldDft<D>& in, RField<D>& out) const
   {
      fft().inverseTransformSafe(in, out);
   }

   /*
   * Apply forward FFT to an array of r-grid fields.
   */
   template <int D>
   void
   FieldIo<D>::convertRGridToKGrid(DArray< RField<D> > const & in,
                                   DArray< RFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   /*
   * Apply forward FFT to a single r-grid field.
   */
   template <int D>
   void
   FieldIo<D>::convertRGridToKGrid(RField<D> const & in,
                                   RFieldDft<D>& out) const
   {
      fft().forwardTransform(in, out);
   }

   /*
   * Test if an RField<D> has declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RField<D> const & in, double epsilon,
                                bool verbose) const
   {
      checkWorkDft();
      fft().forwardTransform(in, workDft_);
      return hasSymmetry(workDft_, epsilon, verbose);
   }

   /*
   * Test if an RFieldDft has the declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RFieldDft<D> const & in, double epsilon,
                                bool verbose) const
   {
      UTIL_CHECK(basis().isInitialized());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> waveCoeff;         // coefficient from wave
      std::complex<double> rootCoeff;         // coefficient from root
      std::complex<double> diff;              // coefficient difference
      int is;                                 // star index
      int iw;                                 // wave index
      int beginId, endId;                     // star begin, end ids
      int rank;                               // dft grid rank of wave

      double cancelledError(0.0);   // max error from cancelled stars
      double uncancelledError(0.0); // max error from uncancelled stars

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(in.dftDimensions());

      // Loop over all stars
      for (is = 0; is < basis().nStar(); ++is) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {

            // Check that coefficients are zero for all waves in star
            beginId = starPtr->beginId;
            endId = starPtr->endId;
            for (iw = beginId; iw < endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!wavePtr->implicit) {
                  rank = dftMesh.rank(wavePtr->indicesDft);
                  waveCoeff = std::complex<double>(in[rank][0], in[rank][1]);
                  if (std::abs(waveCoeff) > cancelledError) {
                     cancelledError = std::abs(waveCoeff);
                     if ((!verbose) && (cancelledError > epsilon)) {
                        return false;
                     }
                  }
               }
            }

         } else {

            // Check consistency of coeff values from all waves
            bool hasRoot = false;
            beginId = starPtr->beginId;
            endId = starPtr->endId;
            for (iw = beginId; iw < endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  rank = dftMesh.rank(wavePtr->indicesDft);
                  waveCoeff = std::complex<double>(in[rank][0], in[rank][1]);
                  waveCoeff /= wavePtr->coeff;
                  if (hasRoot) {
                     diff = waveCoeff - rootCoeff;
                     if (std::abs(diff) > uncancelledError) {
                        uncancelledError = std::abs(diff);
                        if ((!verbose) && (uncancelledError > epsilon)) {
                           return false;
                        }
                     }
                  } else {
                     rootCoeff = waveCoeff;
                     hasRoot = true;
                  }
               }
            }

         }

      } //  end loop over star index is

      if ((cancelledError < epsilon) && (uncancelledError < epsilon)) {
         return true;
      } else if (verbose) {
         Log::file() << std::endl
                     << "Maximum coefficient of a cancelled star: "
                     << cancelledError << std::endl
                     << "Maximum error of coefficient for uncancelled star: "
                     << uncancelledError << std::endl;
      }
      return false;
   }

   // Grid Manipulation Utilties

   template <int D>
   void
   FieldIo<D>::expandFieldsDimension(std::ostream &out,
                                     DArray<RField<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     int d,
                                     DArray<int> newGridDimensions) const

   {
      // Check validity of expanded dimension d and newGridDimensions
      UTIL_CHECK(d > D);
      UTIL_CHECK(d <= 3);
      UTIL_CHECK(newGridDimensions.capacity() == (d - D));

      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Obtain initial dimension of fields
      IntVec<D> meshDimensions = fields[0].meshDimensions();

      int v1 = 1;
      int v2 = 0;
      std::string gname = "";
      FSArray<double, 6> cellParameters;
      if (D == 1){
         //1D expanded to 2D
         if (d == 2){
            DArray<RField<2> > outFields;
            UnitCell<2> cell;
            IntVec<2> dimensions;
            dimensions[0] = meshDimensions[0];
            dimensions[1] = newGridDimensions[0];
            if (dimensions[0] == dimensions[1]){
               cell.set(UnitCell<2>::Square, unitCell.parameters());
            } else {
               cellParameters[0] = unitCell.parameters()[0];
               cellParameters[1] = (double)dimensions[1]/dimensions[0] * cellParameters[0];
               cell.set(UnitCell<2>::Rectangular, cellParameters);
            }

            // Allocate outFields
            outFields.allocate(nMonomer);
            for (int i = 0; i < nMonomer; ++i){
               outFields[i].allocate(dimensions);
               int rank = 0;
               for (int j = 0; j < dimensions[1]; ++j){
                  for (int k = 0; k < dimensions[0]; ++k){
                     outFields[i][rank] = fields[i][k];
                     rank++;
                  }
               }
            }

            // Write Header
            Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gname, nMonomer);
            out << "mesh " <<  std::endl
                << "           " << dimensions << std::endl;

            // Write fields
            MeshIterator<2> itr(dimensions);
            for (itr.begin(); !itr.atEnd(); ++itr) {
               for (int j = 0; j < nMonomer; ++j) {
                  out << "  " << Dbl(outFields[j][itr.rank()], 18, 15);
               }
               out << std::endl;
            }
         }
         
         //1D expanded to 2D
         if (d == 3){
            // Set parameter
            DArray<RField<3> > outFields;
            UnitCell<3> cell;
            IntVec<3> dimensions;
            int rank = 0;
            // Set dimensions
            dimensions[0] = meshDimensions[0];
            dimensions[1] = newGridDimensions[0];
            dimensions[2] = newGridDimensions[1];
            cell.set(UnitCell<3>::Cubic, unitCell.parameters());
            // Assign unit cell
            if (dimensions[2] == dimensions[1]){
               if (dimensions[1] == dimensions[0] ){
                  cell.set(UnitCell<3>::Cubic, unitCell.parameters());
               } else {
                  cellParameters[0] = unitCell.parameters()[0];
                  cellParameters[1] = (double)dimensions[2]/dimensions[0] * cellParameters[0];
               }
            }
            else {
               if (dimensions[1] == dimensions[0]){
                  cellParameters[0] = unitCell.parameters()[0];
                  cellParameters[1] = (double)dimensions[2]/dimensions[0] * cellParameters[0];
                  cell.set(UnitCell<3>::Tetragonal, cellParameters);
               }else if (dimensions[2] == dimensions[0]){
                  cellParameters[0] = unitCell.parameters()[0];
                  cellParameters[1] = (double)dimensions[1]/dimensions[0] * cellParameters[0];
                  cell.set(UnitCell<3>::Tetragonal, cellParameters);
               } else{
                  cellParameters[0] = unitCell.parameters()[0];
                  cellParameters[1] = (double)dimensions[1]/dimensions[0] * cellParameters[0];
                  cellParameters[2] = (double)dimensions[2]/dimensions[0] * cellParameters[0];
                  cell.set(UnitCell<3>::Orthorhombic, cellParameters);
               }

            }
            
            // Allocate outFields
            outFields.allocate(nMonomer);
            for (int i = 0; i < nMonomer; ++i){
               outFields[i].allocate(dimensions);
            }
            for (int l = 0; l < dimensions[2]; ++l){
               for (int k = 0; k< dimensions[1]; ++k){
                  for(int j = 0; j< dimensions[0]; ++j){
                     for (int i = 0; i < nMonomer; ++i){
                        outFields[i][rank] = fields[i][j];
                     }
                     ++rank;
                  }
               }
            }
            
            // Write Header
            Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gname, nMonomer);
            out << "mesh " <<  std::endl
                << "           " << dimensions << std::endl;
                
            // Write fields
            MeshIterator<3> itr(dimensions);
            for (itr.begin(); !itr.atEnd(); ++itr) {
               for (int j = 0; j < nMonomer; ++j) {
                  out << "  " << Dbl(outFields[j][itr.rank()], 18, 15);
               }
               out << std::endl;
            }
         }
      }
      
      //2D expanded to 3D
      if (D == 2){
         // Set parameter
         DArray<RField<3> > outFields;
         UnitCell<3> cell;
         IntVec<3> dimensions;
         int rank = 0;
         // Set dimensions
         dimensions[0] = meshDimensions[0];
         dimensions[1] = meshDimensions[1];
         dimensions[2] = newGridDimensions[0];
         // Set unit cell
         if (unitCell.lattice() == UnitCell<2>::Square){
            if (newGridDimensions[0] == meshDimensions[0]){
               cell.set(UnitCell<3>::Cubic, unitCell.parameters());
            } else {
               cellParameters[0] = unitCell.parameters()[0];
               cellParameters[1] = (double)dimensions[2]/dimensions[0] * cellParameters[0];
               cell.set(UnitCell<3>::Tetragonal, cellParameters);
            }
         }
         if (unitCell.lattice() == UnitCell<2>::Rectangular){
            if (newGridDimensions[0] == meshDimensions[0] || newGridDimensions[0] == meshDimensions[1]){
               cell.set(UnitCell<3>::Tetragonal, unitCell.parameters());
            } else {
               cellParameters[0] = unitCell.parameters()[0];
               cellParameters[1] = unitCell.parameters()[1];
               cellParameters[2] = (double)dimensions[2]/dimensions[0] * cellParameters[0];
               cell.set(UnitCell<3>::Orthorhombic, cellParameters);
            }
         }

         if (unitCell.lattice() == UnitCell<2>::Hexagonal){
            cellParameters[0] = unitCell.parameters()[0];
            cellParameters[1] = (double)dimensions[2]/dimensions[0] * cellParameters[0];
            cell.set(UnitCell<3>::Hexagonal, cellParameters);
         }

         // Allocate outFields
         outFields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i){
            outFields[i].allocate(dimensions);
         }
         
         int q = 0;
         int r = 0;
         int s = 0;
         int n1 =0;
         int n2 =0;
         int n3 =0;
         while (n3 < dimensions[2]) {
            q = 0;
            n2 = 0;
            while (n2 < dimensions[1]) {
               r = q;
               n1 = 0;
               while (n1 < dimensions[0]) {
                  for (int i = 0; i < nMonomer; ++i) {
                     outFields[i][s] = fields[i][r];
                  }
                  r = r + dimensions[1];
                  ++s;
                  ++n1;
               }
               ++q;
               ++n2;
            }
            ++n3;
         }
         // Write Header
         Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gname, nMonomer);
         out << "mesh " <<  std::endl
             << "           " << dimensions << std::endl;

         // Write fields
         MeshIterator<3> itr(dimensions);
         for (itr.begin(); !itr.atEnd(); ++itr) {
            for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(outFields[j][itr.rank()], 18, 15);
            }
            out << std::endl;
         }
      }
   }

   template <int D>
   void
   FieldIo<D>::expandFieldsDimension(std::string filename,
                                     DArray<RField<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     int d,
                                     DArray<int> newGridDimensions) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      expandFieldsDimension(file, fields, unitCell, d, newGridDimensions);
      file.close();
   }

   template <int D>
   void FieldIo<D>::replicateUnitCell(std::ostream &out,
                                      DArray<RField<D> > const & fields,
                                      UnitCell<D> const & unitCell,
                                      IntVec<D> const & replicas) const

   {
      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Obtain initial dimensions of fields
      IntVec<D> meshDimensions = fields[0].meshDimensions();

      // Define dimension of replicated fields
      IntVec<D> replicateDimensions;
      for (int i = 0; i < D; ++i){
         UTIL_CHECK(replicas[i] != 0);
         replicateDimensions[i] = replicas[i] * meshDimensions[i];
      }

      // Set up new Unit Cell
      UnitCell<D> cell;
      FSArray<double, 6> parameters;
      int nParameter = unitCell.nParameter();
      for (int i = 0; i < nParameter; i++){
         parameters[i]=  replicas[i]* unitCell.parameters()[i];
      }
      cell.set(unitCell.lattice(), parameters);

      DArray<RField<D> > outFields;
      // Allocate outFields
      outFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i){
         outFields[i].allocate(replicateDimensions);
      }

      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      int ybeginPtr = 0;
      int zbeginPtr = 0;
      int rank = 0;

      if (D == 1){
         parameters[0]= replicas[0] * unitCell.parameters()[0];
         cell.set(unitCell.lattice(), parameters);
         for (int counter = 0; counter< replicas[0]; counter++){
            for (int j = 0; j < meshDimensions[0]; j++) {
               for (int i = 0; i < nMonomer; i++){
                  outFields[i][rank] = fields[i][j];
               }
               rank++;
            }
         }
      }

      if (D == 2){
         for (int yCounter = 0; yCounter < replicas[0]; yCounter++){
            n2 = 0;
            while (n1 < meshDimensions[0]){
               ybeginPtr = n1 * meshDimensions[1];
               for (int xCounter = 0; xCounter< replicas[1]; xCounter++){
                  n1 = 0;
                  int r = ybeginPtr;
                  while (n2 < meshDimensions[1]){
                     for (int i = 0; i < nMonomer; i++){
                        outFields[i][rank] = fields[i][r];
                     }
                     rank++;
                     n2++;
                     r++;
                  }
               }
               n1++;
            }
         }
      }

      if (D == 3){
         for (int zCounter = 0; zCounter < replicas[0]; zCounter++){
            n1 = 0;
            while (n1 < meshDimensions[0]){
               zbeginPtr =  n1* meshDimensions[2] * meshDimensions[1];
               for (int yCounter = 0; yCounter < replicas[1]; yCounter++){
                  n2 = 0;
                  while (n2 < meshDimensions[1]){
                     ybeginPtr = zbeginPtr + n2 * meshDimensions[2];
                     for (int xCounter = 0; xCounter< replicas[2]; xCounter++){
                        n3 = 0;
                        int r = ybeginPtr;
                        while (n3 < meshDimensions[2]){
                           for (int i = 0; i < nMonomer; i++){
                              outFields[i][rank] = fields[i][r];
                           }
                           rank++;
                           n3++;
                           r++;
                        }
                     }
                     n2++;
                  }
               }
               n1++;
            }
         }
      }

      // Write Header
      int v1 = 1;
      int v2 = 0;
      std::string gname = "";
      Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gname, nMonomer);
      out << "mesh " <<  std::endl
          << "           " << replicateDimensions << std::endl;

      // Write Fields
      FieldIo<D>::writeFieldsRGrid(out, outFields, cell, false, false, false);
   }

   template <int D>
   void FieldIo<D>::replicateUnitCell(std::string filename,
                                      DArray<RField<D> > const & fields,
                                      UnitCell<D> const & unitCell,
                                      IntVec<D> const & replicas) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      replicateUnitCell(file, fields, unitCell, replicas);
      file.close();
   }

   // Private utility functions

   template <int D>
   void FieldIo<D>::checkWorkDft() const
   {
      if (!workDft_.isAllocated()) {
         workDft_.allocate(mesh().dimensions());
      } else {
         UTIL_CHECK(workDft_.meshDimensions() == fft().meshDimensions());
      }
   }

} // namespace Rpc
} // namespace Pscf
#endif

#ifndef PRDC_CL_FIELD_IO_H
#define PRDC_CL_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>     // nested LatticeSystem enum
#include <pscf/math/IntVec.h>          // template with default
#include <util/containers/DArray.h>    // member

// Forward declarations 
namespace Util {
   class FileMaster;
}
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Prdc {
namespace Cl {

   using namespace Util;
   using namespace Pscf;
   using namespace Pscf::Prdc;

   /**
   * File input/output, format conversions and other utilities for fields.
   *
   * This class template provides functions to read and write real-valued
   * fields in any of three representations (symmetry-adapted basis,
   * real-space r-grid, or Fourier space k-grid form), functions to
   * convert among these three representations, and other utilities for
   * manipulating fields and field files.
   *
   * <b>Template parameters:</b>
   *
   *    - D    : integer dimension of space, i.e., 1, 2, or 3
   *    - CFT  : complex field type, e.g., CField<D>
   *    - FFT  : fast Fourier transform type, e.g., FFT<D>
   *
   * <b>Subclasses:</b>
   * The Cl:FieldIo template is a base class for two class templates 
   * named FieldIo that are defined in namespaces Pscf::Cpc and Pscf::Cpg.
   * The Pscf::Cpc::FieldIo<int D> template is derived from a partial
   * specialization of Cl::FieldIo with parameters CFT = Cpu::CField<D> 
   * and FFT = Cpu::FFT<D> that are defined in the Prdc::Cpu namespace,
   * and that use standard CPU hardware.  The analogous class template 
   * Cpg::FieldIo<int D> in the Pscf::Cpg namespace is derived from a 
   * partial specialization of Cl::FieldIo in which these two parameters 
   * are class templates with the same names (CField and FFT) that are 
   * defined in the Prdc::Cuda namespace and use GPU hardware.
   *
   * <b> Pure virtual member functions </b>: This class template defines
   * several pure virtual functions for which different implementations
   * are required for Cpu and Cuda code. Cpu and Cuda implementations
   * of these functions, which are defined in the Cpc::FieldIo<D> and
   * Cpg::FieldIo<D> subclasses, differ because the Cuda versions must
   * explicitly transfer data between Cpu and Gpu memory.
   *
   * \ingroup Prdc_Cl_Module
   */
   template <int D, class CFT, class FFT>
   class FieldIo
   {

   public:

      /// \name Construction, Initialization and Destruction
      ///@{

      /**
      * Constructor.
      */
      FieldIo();

      /**
      * Destructor.
      */
      virtual ~FieldIo();

      /**
      * Create associations with other members of the parent Domain.
      *
      * This function may be called within the constructor of the
      * Domain object, since addresses of other members of the Domain
      * are known at this point.
      *
      * \param mesh  associated spatial discretization Mesh<D>
      * \param fft   associated FFT object for fast transforms
      * \param lattice  lattice system type (enumeration value)
      */
      void associate(Mesh<D> const & mesh,
                     FFT const & fft,
                     typename UnitCell<D>::LatticeSystem const & lattice);

      /**
      * Create an association with a FileMaster.
      *
      * The FileMaster is used to open and close files in all member
      * functions that take file name arguments and that open and close
      * files. This allows prefixes for input and output files (if any)
      * to be automatically prepended to file names.
      *
      * \param fileMaster  associated FileMaster (for file paths)
      */
      void setFileMaster(FileMaster const & fileMaster);

      /**
      * Set the number of monomer types.
      *
      * This is used to allocate private arrays of fields used by some
      * functions that read field data into this workspace memory. This
      * function may only be called once, shortly after the value of
      * nMonomer is read from the parameter file.
      *
      * \param nMonomer  number of monomer types
      */
      void setNMonomer(int nMonomer);

      ///@}
      /// \name Field File IO - Real Space Grid Format
      ///@{

      /**
      * Read array of complex r-grid fields from an input stream.
      *
      * Upon successful return, element fields[i] of the fields array is
      * the instance of CFT containing the r-grid field associated with
      * monomer type i.
      *
      * On entry, array fields must either be unallocated or be allocated
      * with capacity equal to the number of monomers in the field file,
      * with mesh dimensions for each field equal to mesh().dimensions().
      * If it is unallocated, it will be allocated within this function.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of r-grid fields
      * \param unitCell  associated crystallographic unit cell
      */
      virtual
      void readFields(std::istream& in,
                      DArray<CFT>& fields,
                      UnitCell<D> & unitCell) const = 0;

      /**
      * Read an array of r-grid fields from a named file.
      *
      * This function opens an input file with the specified filename,
      * reads fields in real-space grid format from that file, and
      * then closes the file. The overloaded readFieldsRGrid member
      * function that takes a std::istream& argument is called to read
      * the file.
      *
      * \param filename  name of input file
      * \param fields  array of r-grid fields (instances of CFT)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFields(std::string filename,
                      DArray<CFT>& fields,
                      UnitCell<D> & unitCell) const;

      /**
      * Read data for array of r-grid fields, with no header section.
      *
      * This function reads the data section of the rgrid-field format,
      * without first reading a header.
      *
      * On entry, array fields must either be unallocated or be allocated
      * with capacity equal to the number of monomers in the field file,
      * with mesh dimensions for each field equal to mesh().dimensions().
      * If it is unallocated, it will be allocated within this function.
      *
      * \param in  input file stream
      * \param fields  array of r-grid fields (instances of CFT)
      * \param nMonomer  expected number of monomer types (input
      */
      virtual
      void readFieldsdData(std::istream& in,
                           DArray<CFT>& fields,
                           int nMonomer) const = 0;

      /**
      * Read a single r-grid field from an input stream.
      *
      * On entry, the field must either be unallocated or be allocated
      * with a mesh dimension equal to mesh().dimensions(). If it is
      * unallocated, it will be allocated within this function.
      *
      * \param in  input stream (i.e., input file)
      * \param field  single r-grid field (instance of CFT)
      * \param unitCell  associated crystallographic unit cell
      */
      virtual
      void readField(std::istream &in,
                     CFT & field,
                     UnitCell<D>& unitCell) const = 0;

      /**
      * Read a single field from a named file.
      *
      * This function opens an input file with the specified filename,
      * reads a field in grid format from that file, and then closes
      * the file. The overloaded readField member function that takes 
      * a std::istream& argument is called to read the file.
      *
      * \param filename  name of input file
      * \param field  single field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      */
      void readField(std::string filename,
                     CFT & field,
                     UnitCell<D>& unitCell) const;

      /**
      * Write an array of fields to an output stream.
      *
      * On entry, the container fields must be allocated, and the mesh
      * dimensions of each field must equal mesh().dimensions(). The
      * writeHeader argument may be set false to completely suppress
      * writing of the file header. The isSymmetric argument is only
      * relevant if writeHeader is true.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  iff true, write file header
      * \param writeMeshSize Should mesh size be written at end of header?
      */
      virtual
      void writeFields(std::ostream& out,
                       DArray<CFT> const & fields,
                       UnitCell<D> const & unitCell,
                       bool writeHeader = true,
                       bool writeMeshSize = true) const = 0;

      /**
      * Write an array of r-grid fields to a named file.
      *
      * This function opens a file, writes field file header and data
      * to the file, and closes the file. The overloaded writeFieldsGrid
      * member function that takes a std::ostream& argument is called
      * internally with writeHeader == true to  write the data.
      *
      * \param filename  name of output file
      * \param fields  array of CFT objects (fields on a grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFields(std::string filename,
                       DArray<CFT> const & fields,
                       UnitCell<D> const & unitCell) const;

      /**
      * Write a single r-grid field to an an output stream.
      *
      * On entry, the field container must be allocated with mesh
      * dimensions equal to mesh().dimensions(). The writeHeader
      * argument may be set false to suppress writing of the file
      * header.
      *
      * \param out  output stream
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  should a file header be written?
      */
      virtual
      void writeField(std::ostream &out,
                      CFT const & field,
                      UnitCell<D> const & unitCell,
                      bool writeHeader = true) const = 0;

      /**
      * Write a single r-grid field to a named file.
      *
      * This function opens a file, writes the header and data for a
      * single field to the file, and closes that file. The overloaded
      * writeFieldGrid member function that takes a std::ostream&
      * argument is called internally with writeHeader == true to
      * write the file.
      *
      * \param filename  name of output file
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      */
      void writeField(std::string filename,
                      CFT const & field,
                      UnitCell<D> const & unitCell) const;

      ///@}
     
      #if 0
      /// \name Field Array Format Conversion
      ///@{

      /**
      * Convert a single field from k-grid to r-grid format.
      *
      * This function simply calls the inverse FFT for a single field.
      *
      * \param in  field in discrete Fourier format (k-grid)
      * \param out  field defined on real-space grid (r-grid)
      */
      void convertKGridToRGrid(CFT const & in,
                               CFT & out) const;

      /**
      * Convert an array of fields from r-grid to k-grid (Fourier) format.
      *
      * This function simply calls the forward FFT repeatedly for an
      * array of fields.
      *
      * \param in  fields defined on real-space grid (r-grid)
      * \param out  fields in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(DArray<CFT> const & in,
                               DArray<CFT> & out) const;

      /**
      * Convert a field from r-grid to k-grid (Fourier) format.
      *
      * This function simply calls the forward FFT for a single field.
      *
      * \param in   field defined on real-space grid (r-grid)
      * \param out  field in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(CFT const & in, CFT & out) const;

      ///@}
      /// \name Field File Format Conversion
      ///@{

      /**
      * Convert a field file from Fourier (k-grid) to r-grid format.
      *
      * The number of monomers in the input file must equal the number
      * set by the setNMonomer(int) member function.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (r-grid format)
      */
      void convertKGridToRGrid(std::string const & inFileName,
                               std::string const & outFileName) const;

      /**
      * Convert a field file from r-grid to Fourier (k-grid) format.
      *
      * The number of monomers in the input file must equal the number
      * set by the setNMonomer(int) member function.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (k-grid format)
      */
      void convertRGridToKGrid(std::string const & inFileName,
                               std::string const & outFileName) const;

      ///@}
      /// \name Test Field Equality
      ///@{

      /**
      * Compare two fields in r-grid form, write a report to Log file.
      *
      * Outputs maximum and root-mean-squared differences to the
      * standard Log file.
      *
      * \param field1  first array of fields (r-grid form)
      * \param field2  second array of fields (r-grid form)
      */
      virtual
      void compareFields(DArray< CFT > const & field1,
                         DArray< CFT > const & field2) const = 0;

      /**
      * Compare two r-grid field files, write a report to Log file.
      *
      * \param filename1  name of 1st field file
      * \param filename2  name of 2nd field file
      */
      void compareFields(std::string const & filename1,
                         std::string const & filename2) const;

      ///@}
      #endif

      /// \name Field File IO Utilities
      ///@{

      /**
      * Reader header of field file (fortran PSCF format)
      *
      * This reads the common part of the header for all PSCF field file
      * formats. This contains the dimension of space, the lattice
      * system, a list of unit cell parameters, the space group name as
      * an optional parameter, and the number of monomer types. The unit
      * cell data is used to update a UnitCell<D> that is passed as a
      * parameter.
      *
      * On return, parameter nMonomer contains the number of monomer
      * types declared in the field file header.  This function does
      * not require the number of monomers declared in the field file
      * header to match the value of nMonomer set by the setNMonomer
      * member function.
      *
      * Consistency checks (Exceptions thrown on failure):
      *
      * The value of "dim" in the header file must match the template
      * parameter D.  If the UnitCell<D> object passed to this function
      * already contains a non-null lattice type, it must match the
      * lattice system in the header file.
      *
      * \param in  input stream (i.e., file)
      * \param nMonomer  number of monomer types in the header (output)
      * \param unitCell  associated crystallographic unit cell (output)
      */
      void readFieldHeader(std::istream& in,
                           int& nMonomer,
                           UnitCell<D> & unitCell) const;

      /**
      * Write header for field file (fortran pscf format).
      *
      * \param out  output stream (i.e., file)
      * \param nMonomer  number of monomer types or fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Should a space group be declared?
      */
      void writeFieldHeader(std::ostream& out,
                            int nMonomer,
                            UnitCell<D> const & unitCell) const;

      ///@}
      /// \name Accessor functions (const references)
      ///@{

      /**
      * Get spatial discretization mesh by const reference.
      */
      Mesh<D> const & mesh() const
      {
         UTIL_ASSERT(meshPtr_);
         return *meshPtr_;
      }

      /**
      * Get associated FileMaster by const reference.
      */
      FileMaster const & fileMaster() const
      {
         UTIL_ASSERT(fileMasterPtr_);
         return *fileMasterPtr_;
      }

      ///@}

   protected:

      /**
      * Get the lattice type enum value by const reference.
      */
      typename UnitCell<D>::LatticeSystem const & lattice() const
      {
         UTIL_ASSERT(latticePtr_);
         return *latticePtr_;
      }

      /**
      * Get FFT object by const reference
      */
      FFT const & fft() const
      {
         UTIL_ASSERT(fftPtr_);
         return *fftPtr_;
      }

      /**
      * Check if r-grid workspace is allocated, allocate if necessary.
      */
      void checkAllocateGrid() const;

   private:

      /// Pointer to spatial discretization mesh
      Mesh<D> const * meshPtr_;

      /// Pointer to FFT object
      FFT const * fftPtr_;

      /// Pointer to lattice system
      typename UnitCell<D>::LatticeSystem const * latticePtr_;

      /// Pointer to Filemaster (holds paths to associated I/O files)
      FileMaster const * fileMasterPtr_;

      /// Number of monomer types
      int nMonomer_;

      // Mutable work space used by functions that read field files

      /// Work array of complex fields on a mesh.
      mutable DArray<CFT> tmpFields_;

      /// Is tmpFields_ allocated?
      mutable bool isAllocated_;

   };

} // namespace Cl
} // namespace Prdc
} // namespace Pscf
#endif

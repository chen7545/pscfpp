#ifndef CPC_FIELD_IO_H
#define CPC_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cl/FieldIo.h>      // base class template
#include <prdc/cpu/CField.h>      // template parameter
#include <prdc/cpu/FFT.h>         // template parameter

// Forward declarations
namespace Util {
   class FileMaster;
   template <typename T> class DArray;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

// Explicit instantiation declaration for base class
namespace Pscf {
   namespace Prdc {
      extern template class Cl::FieldIo<1, Cpu::CField<1>, Cpu::FFT<1> >;
      extern template class Cl::FieldIo<2, Cpu::CField<2>, Cpu::FFT<2> >;
      extern template class Cl::FieldIo<3, Cpu::CField<3>, Cpu::FFT<3> >;
   }
}

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * File input/output operations for fields.
   *
   * Please refer to the documentation of the base class template
   * Prdc::Cl::FieldIo for complete API documentation. The public
   * interface of this class is identical to that of the base class.
   *
   * This class template is derived from a partial specialization of
   * the template Prdc::Cl::FieldIo<D, CFT, FFT> using classes
   * CFT = CField<D> and FFT = FFT<D> that are all defined in the
   * Prdc::Cpu subspace, and that all use conventional CPU hardware.
   * An analogous class template named Rpg::FieldIo that is defined
   * in the Pscf::Rpg namespace instead uses a GPU.
   *
   * The member functions defined in this class are all implementations of
   * pure virtual functions declared in the base class, Prdc::Cl::FieldIo.
   * These are all functions for which different implementations are
   * required for the CPU and GPU variants, usually because the GPU
   * implementation requires data transfer between host and device.
   *
   * \ingroup Cpc_Field_Module
   */
   template <int D>
   class FieldIo
     : public  Cl::FieldIo< D, CField<D>, FFT<D> >
   {

   public:

      /**
      * Read array of CField objects (r-grid fields) from a stream.
      *
      * See documentation of analogous function in Prdc::Cl::FieldIo.
      *
      * \param in  input file stream
      * \param fields  array of complex fields
      * \param unitCell  associated crystallographic unit cell
      */
      void readFields(std::istream& in,
                      DArray< CField<D> >& fields,
                      UnitCell<D> & unitCell) const override;

      /**
      * Read data for an array of r-grid fields, with no header section.
      *
      * See documentation of analogous function in Prdc::Cl::FieldIo.
      *
      * \param in  input file stream
      * \param fields  array of complex fields 
      * \param nMonomer  number of monomer types
      */
      void readFieldsData(std::istream& in,
                          DArray< CField<D> >& fields,
                          int nMonomer) const override;

      /**
      * Read a single CField (field on an r-space grid) from a stream.
      *
      * See documentation of analogous function in Prdc::Cl::FieldIo.
      *
      * \param in  input file stream
      * \param field  fields defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      */
      void readField(std::istream &in,
                     CField<D> & field,
                     UnitCell<D>& unitCell) const override;

      /**
      * Write array of CField objects (fields on r-space grid) to a stream.
      *
      * See documentation of analogous function in Prdc::Cl::FieldIo.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of CField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write file header if true
      * \param writeMeshSize  Should mesh size be written in header?
      */
      void writeFields(std::ostream& out,
                       DArray< CField<D> > const & fields,
                       UnitCell<D> const & unitCell,
                       bool writeHeader = true,
                       bool writeMeshSize = true) const override;

      /**
      * Write a single CField (field on an r-space grid) to a stream.
      *
      * See documentation of analogous function in Prdc::Cl::FieldIo.
      *
      * \param out  output stream
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  Should a file header be written?
      */
      void writeField(std::ostream &out,
                      CField<D> const & field,
                      UnitCell<D> const & unitCell,
                      bool writeHeader = true)
      const override;

      #if 0
      /**
      * Compare two arrays of complex fields, output a report.
      *
      * Outputs maximum and root-mean-squared differences to the standard
      * Log file.
      *
      * \param field1  first array of fields (r-grid format)
      * \param field2  second array of fields (r-grid format)
      */
      void compareFields(DArray< CField<D> > const & field1,
                         DArray< CField<D> > const & field2)
      const override;
      #endif

      /// Alias for base class
      using Base = Prdc::Cl::FieldIo< D, CField<D>, FFT<D> >;

      // Inherited public member functions
      using Base::associate;
      using Base::setFileMaster;
      using Base::readFields;
      using Base::readFieldsData;
      using Base::readField;
      using Base::writeFields;
      using Base::writeField;
      //using Base::convertKGridToRGrid;
      //using Base::convertRGridToKGrid;
      //using Base::compareFieldsRGrid;
      using Base::readFieldHeader;
      using Base::writeFieldHeader;
      using Base::mesh;
      using Base::fileMaster;

   protected:

      // Inherited protected member functions
      using Base::lattice;
      using Base::fft;

   };

   extern template class FieldIo<1>;
   extern template class FieldIo<2>;
   extern template class FieldIo<3>;

} // namespace Cpc

} // namespace Pscf
#endif

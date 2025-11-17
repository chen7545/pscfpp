#ifndef PRDC_CFIELD_IO_UTIL_H
#define PRDC_CFIELD_IO_UTIL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>   // Template with a default parameter

// Forward class declarations 
namespace Util {
   template <typename T> class DArray;
}
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
      template <int D> class Basis;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf;

   // Templates for complex field data IO

   /**
   * Read data for an array of complex fields, with no header section.
   *
   * This function reads the data section of a complex field file for
   * multiple monomer types, with no header section.
   *
   * The template parameter AT must be an array type that provides an
   * overloaded [] subscript operator that returns a complex number of
   * some type CT. The function template <typename CT, typename RT> 
   * assign(CT&, RT const&, RT const&) must be defined. Valid complex
   * array types include DArray<fftw_complex> and HostDArray<cudaReal>.
   *
   * \ingroup Prdc_Cl_Module
   *
   * \param in  input file stream
   * \param fields  array of complex fields (out)
   * \param nMonomer  number of monomer types (in)
   * \param dimensions  vector of mesh dimensions (in)
   */
   template <int D, class AT>
   void readCFieldData(std::istream& in,
                      DArray< AT >& fields,
                      int nMonomer,
                      IntVec<D> const& dimensions);
 
   /**
   * Read data for a single r-grid field, with no header section.
   *
   * This function reads the data section of a complex field file for
   * a single monomer type, with no header section.
   *
   * The template parameter AT must be an array type that provides an
   * overloaded [] subscript operator that returns a complex number.
   *
   * \ingroup Prdc_Cl_Module
   *
   * \param in  input file stream
   * \param field  array containing a single complex field (out)
   * \param dimensions  vector of mesh dimensions (in)
   */
   template <int D, class AT>
   void readCFieldData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dimensions);
 
   /**
   * Write data for array of r-grid fields, with no header section.
   *
   * This function writes the data section of the complex field file
   * format for a multiple monomer types, with no header section.
   *
   * The template parameter AT must be an array type that provides 
   * an overloaded [] subscript operator that returns a complex number.
   *
   * \ingroup Prdc_Cl_Module
   *
   * \param out  output file stream
   * \param fields  array of complex fields (in)
   * \param nMonomer  number of monomer types (in)
   * \param dimensions  vector of mesh dimensions (in)
   */
   template <int D, class AT>
   void writeCFieldData(std::ostream& out, 
                       DArray<AT> const& fields,
                       int nMonomer,
                       IntVec<D> const& dimensions);
 
   /**
   * Write data for a single r-grid field, with no header section.
   *
   * This function writes the data section of complex field file format
   * for a single monomer type, with no header.
   *
   * The template parameter AT must be an array type that provides 
   * an overloaded [] subscript operator that returns a complex number.
   *
   * \ingroup Prdc_Cl_Module
   *
   * \param out  output file stream
   * \param field  array containing a single complex field (out)
   * \param dimensions  vector of mesh dimensions
   */
   template <int D, class AT>
   void writeCFieldData(std::ostream& out, 
                       AT const& field,
                       IntVec<D> const& dimensions);

} // namespace Prdc
} // namespace Pscf
#endif

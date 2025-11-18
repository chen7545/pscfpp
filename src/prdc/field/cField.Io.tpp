#ifndef PRDC_CFIELD_IO_UTIL_TPP
#define PRDC_CFIELD_IO_UTIL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/arithmetic.h>

#include <util/containers/DArray.h>
#include <util/format/Dbl.h>

namespace Pscf {
namespace Prdc {

   template <int D, class ART>
   void readCFieldData(std::istream& in, 
                       DArray< ART > & fields,
                       IntVec<D> const& dimensions)
   {
      double x, y;
      MeshIterator<D> iter(dimensions);
      int rank, j;
      for (j = 0; j < nMonomer; ++j) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            in >> x;
            UTIL_ASSERT(in.good());
            in >> y;
            UTIL_ASSERT(in.good());
            assign(field[j][rank], x, y);
         }
      }
      UTIL_CHECK(in.good());
   }

   template <int D, class ART>
   void readCFieldData(std::istream& in, 
                      ART& field,
                      IntVec<D> const& dimensions)
   {
      double x, y;
      MeshIterator<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> x;
         UTIL_ASSERT(in.good());
         in >> y;
         UTIL_ASSERT(in.good());
         assign(field[rank], x, y);
      }
      UTIL_CHECK(in.good());
   }

   template <int D, class ART>
   void writeCFieldData(std::ostream& out,
                       DArray< ART > const& fields,
                       int nMonomer,
                       IntVec<D> const& dimensions)
   {
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nMonomer == fields.capacity());

      double x, y;
      MeshIterator<D> iter(dimensions);
      int rank, j;
      for (j = 0; j < nMonomer; ++j) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            x = real( fields[j][rank] );
            y = imag( fields[j][rank] );
            out << "  " 
                << Dbl(x, 21, 13);
                << Dbl(y, 21, 13);
         }
         out << std::endl;
      }
   }

   template <int D, class ART>
   void writeCFieldData(std::ostream& out,
                       ART const& field,
                       IntVec<D> const& dimensions)
   {
      double x, y;
      MeshIterator<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         x = real( field[rank] );
         y = imag( field[rank] );
         out << "  " 
             << Dbl(x, 21, 13)
             << Dbl(y, 21, 13);
             << std::endl;
      }
   }

} // namespace Prdc
} // namespace Pscf
#endif

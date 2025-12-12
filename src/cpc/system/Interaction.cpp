/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Interaction.h"
#include <cmath>

namespace Pscf {
namespace Cpc {

   using namespace Util;

   /*
   * Constructor.
   */
   Interaction::Interaction(bool isCompressible)
    : chi_(),
      zeta_(-1.0),
      nMonomer_(0),
      isCompressible_(isCompressible),
      isInitialized_(false)
   {  setClassName("Interaction"); }

   /*
   * Destructor.
   */
   Interaction::~Interaction()
   {}

   /*
   * Set the number of monomer types.
   */
   void Interaction::setNMonomer(int nMonomer)
   {  
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer; 
      chi_.allocate(nMonomer, nMonomer);
   }

   /*
   * Read parameters from file.
   */
   void Interaction::readParameters(std::istream& in)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(!isInitialized_);
      readDSymmMatrix(in, "chi", chi_, nMonomer());
      if (isCompressible_) {
         //isCompressible_ = readOptional(in, "zeta", zeta_).isActive();
         auto& param = readOptional(in, "zeta", zeta_);
         isCompressible_ = param.isActive();
      }
      isInitialized_ = true;
   }

   /*
   * Set a single Flory-Huggins chi parameter.
   */
   void Interaction::setChi(int i, int j, double chi)
   {
      UTIL_CHECK(isInitialized_);
      chi_(i,j) =  chi; 
      if (i != j) {
         chi_(j,i) = chi;
      }
   }

} // namespace Cpc
} // namespace Pscf

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Interaction.h"
#include <cmath>

namespace Pscf {
namespace Cp {

   using namespace Util;

   /*
   * Constructor.
   */
   Interaction::Interaction()
    : chi_(),
      zeta_(-1.0),
      range_(-1.0),
      alpha_(1.0),
      nMonomer_(0)
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
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer; 
      chi_.allocate(nMonomer, nMonomer);
   }

   /*
   * Read parameters from file.
   */
   void Interaction::readParameters(std::istream& in)
   {
      UTIL_CHECK(nMonomer() > 0);
      readDSymmMatrix(in, "chi", chi_, nMonomer());
      read(in, "zeta", zeta_);
      read(in, "range", range_);

      alpha_ = -0.5/(range_ * range_);
   }

   /*
   * Set a single Flory-Huggins chi parameter.
   */
   void Interaction::setChi(int i, int j, double chi)
   {
      chi_(i,j) =  chi; 
      if (i != j) {
         chi_(j,i) = chi;
      }
   }

   /*
   * Compute & return Fourier-space damping factor for pair interactions.
   */
   double Interaction::g(double kSq) const
   {
      return std::exp(alpha_ * kSq);
   }

} // namespace Cp
} // namespace Pscf

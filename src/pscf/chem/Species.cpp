/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"
#include <cmath>

namespace Pscf { 

   using namespace Util;

   /*
   * Default constructor.
   */
   Species::Species()
    : phi_(0.0),
      mu_(0.0),
      q_(0.0),
      ensemble_(Ensemble::Closed)
   {  setClassName("Species"); }

   /*
   * Read phi or mu (but not both) and set ensemble.
   */
   void Species::readParameters(std::istream& in)
   {
      // Read phi or mu (but not both)
      bool hasPhi = readOptional(in, "phi", phi_).isActive();
      if (hasPhi) {
         ensemble_ = Ensemble::Closed;
      } else {
         ensemble_ = Ensemble::Open;
         read(in, "mu", mu_);
      }
   }

   /*
   *  Set volume fraction (if ensemble is closed).
   */ 
   void Species::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Ensemble::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi;
   }

   /*
   *  Set chemical potential (if ensemble is open).
   */ 
   void Species::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Ensemble::Open);  
      mu_ = mu;
   }

   /*
   * Set q and compute mu or phi (depending on ensemble).
   */
   void Species::setQ(double q)
   {
      q_ = q;
      if (ensemble() == Ensemble::Closed) {
         mu_ = std::log(phi_/q_);
      } else
      if (ensemble() == Ensemble::Open) {
         phi_ = std::exp(mu_)*q_;
      }
   }

} // namespace Pscf

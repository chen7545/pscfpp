/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"

namespace Pscf {
   namespace Prdc { 
      template 
      class Cl::Mixture<1, Cpc::Polymer<1>, Cpc::Solvent<1>, Cpc::Types<1> >;
      template 
      class Cl::Mixture<2, Cpc::Polymer<2>, Cpc::Solvent<2>, Cpc::Types<2> >;
      template 
      class Cl::Mixture<3, Cpc::Polymer<3>, Cpc::Solvent<3>, Cpc::Types<3> >;
   }
   namespace Cpc { 
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}

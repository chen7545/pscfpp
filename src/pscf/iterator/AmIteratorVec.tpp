#ifndef PSCF_AM_ITERATOR_VEC_TPP
#define PSCF_AM_ITERATOR_VEC_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorVec.h"
#include <util/global.h> 

namespace Pscf
{

   using namespace Util;

   /*
   * Vector assignment, a = b .
   */
   template <class Iterator, class T>
   void AmIteratorVec<Iterator,T>::setEqual(T& a, T const & b)
   {  VecOp::assign(a, b); }

   /*
   * Compute and return the inner product of two vectors
   */
   template <class Iterator, class T> 
   double 
   AmIteratorVec<Iterator,T>::dotProduct(T const & a, T const & b)
   {  return VecOp::innerProduct(a, b); }

   /*
   * Compute and return the maximum magnitude element of a vector.
   */
   template <class Iterator, class T>
   double AmIteratorVec<Iterator,T>::maxAbs(T const & a)
   {  return Reduce::maxAbs(a); }

   /*
   * Compute the vector difference a = b - c .
   */
   template <class Iterator, class T>
   void AmIteratorVec<Iterator,T>::subVV(T& a, T const & b, T const & c)
   {  VecOp::subVV(a, b, c); }

   /*
   * Composite a += b*c for vectors a and b, scalar c .
   */
   template <class Iterator, class T>
   void AmIteratorVec<Iterator,T>::addEqVc(T& a, T const & b, double c)
   {  VecOp::addEqVc(a, b, c); }

}
#endif

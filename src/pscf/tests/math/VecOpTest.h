#ifndef PSCF_VEC_OP_TEST_H
#define PSCF_VEC_OP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cpu/VecOp.h>
#include <util/math/Constants.h>
#include <cmath>

using namespace Util;
using namespace Pscf;

class VecOpTest : public UnitTest
{

private:

   // Error tolerance for array equality
   constexpr static double tolerance_ = 1E-10;

   // Array size, large enough to require multiple blocks
   const static int n = 2048; 

   // Input and output arrays, real and complex
   DArray<double> inReal, inReal2, outReal;
   DArray<double> refOutReal;

   // Input scalars, real and complex
   double scalarReal;

   // Input and outputs using standard types, for comparison

public:

   void setUp()
   {
      // Allocate arrays
      inReal.allocate(n);
      inReal2.allocate(n);
      outReal.allocate(n);
      refOutReal.allocate(n);

      // Define "in" arrays with arbitrary data between -1 and 1
      double halfPi = Constants::Pi;
      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;
         inReal[i] = sin(halfPi * frac);
         inReal2[i] = cos(halfPi * frac); 
      }

      // Define input scalars with arbitrary values
      scalarReal = 0.633; 
   }

   void tearDown()
   {}

   // Test VecOp::eqV and VecOp::eqS
   void testEq()
   {
      printMethod(TEST_FUNC);

      VecOp::eqV(outReal, inReal);
      checkEqualReal(outReal, inReal);

      VecOp::eqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   // Test VecOp::addVV and VecOp::addVS
   void testAdd()
   {
      printMethod(TEST_FUNC);

      VecOp::addVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] + inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::addVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] + scalarReal;
      }
      checkEqualReal(outReal, refOutReal);
   }

   void testSub()
   {
      printMethod(TEST_FUNC);

      VecOp::subVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] - inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::subVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] - scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   void testMul()
   {
      printMethod(TEST_FUNC);

      VecOp::mulVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] * inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::mulVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] * scalarReal;
      }
      checkEqualReal(outReal, refOutReal);
   }

   void testDiv()
   {
      printMethod(TEST_FUNC);

      VecOp::divVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] / inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::divVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i] / scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   void testExp()
   {
      printMethod(TEST_FUNC);

      VecOp::expV(outReal, inReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = exp(inReal[i]);
      }
      checkEqualReal(outReal, refOutReal);

   }

   void testAddEq()
   {
      printMethod(TEST_FUNC);

      VecOp::eqV(outReal, inReal);
      VecOp::addEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] += inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::eqV(outReal, inReal);
      VecOp::addEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] += scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   // Test VecOp::subEqV and VecOp::subEqS
   void testSubEq()
   {
      printMethod(TEST_FUNC);

      VecOp::eqV(outReal, inReal);
      VecOp::subEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] -= inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::eqV(outReal, inReal);
      VecOp::subEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] -= scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   void testMulEq()
   {
      printMethod(TEST_FUNC);

      VecOp::eqV(outReal, inReal);
      VecOp::mulEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] *= inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::eqV(outReal, inReal);
      VecOp::mulEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] *= scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   void testDivEq()
   {
      printMethod(TEST_FUNC);

      VecOp::eqV(outReal, inReal);
      VecOp::divEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] /= inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::eqV(outReal, inReal);
      VecOp::divEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) { 
         refOutReal[i] = inReal[i];
         refOutReal[i] /= scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

   }

   // Utility function

   void checkEqualReal(DArray<double>& a, DArray<double>& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a[i] - b[i]) < tolerance_);
      }
   }

};

TEST_BEGIN(VecOpTest)
TEST_ADD(VecOpTest, testEq)
TEST_ADD(VecOpTest, testAdd)
TEST_ADD(VecOpTest, testSub)
TEST_ADD(VecOpTest, testMul)
TEST_ADD(VecOpTest, testDiv)
TEST_ADD(VecOpTest, testExp)
TEST_ADD(VecOpTest, testAddEq)
TEST_ADD(VecOpTest, testSubEq)
TEST_ADD(VecOpTest, testMulEq)
TEST_ADD(VecOpTest, testDivEq)
TEST_END(VecOpTest)

#endif

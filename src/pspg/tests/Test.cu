/*
* This program runs all unit tests in the pspg/tests directory.
*/ 

#include <test/CompositeTestRunner.h>

#include "field/FieldTestComposite.h"
#include "solvers/SolverTestComposite.h"
#include "system/SystemTestComposite.h"
#include "sweep/SweepTestComposite.h"
#include "simulate/SimulatorTest.h"
#include <util/global.h>

TEST_COMPOSITE_BEGIN(PspgNsTestComposite)
addChild(new FieldTestComposite, "field/");
addChild(new SolverTestComposite, "solvers/");
addChild(new SystemTestComposite, "system/");
addChild(new SweepTestComposite, "sweep/");
addChild(new TEST_RUNNER(SimulatorTest), "simulate/");
TEST_COMPOSITE_END

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   try {

      PspgNsTestComposite runner;

      // Add any file prefix given as command line argument
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
      }
   
      // Run all unit test methods
      int failures = runner.run();

      if (failures != 0) {
         failures = 1;
      }
      return failures;

   } catch (...) {

      std::cerr << "Uncaught exception in pspg/tests/Test.cc" << std::endl;
      return 1;

   }
}

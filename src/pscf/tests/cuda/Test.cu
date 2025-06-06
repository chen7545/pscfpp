/*
* This program runs all unit tests in the pscf/tests/math directory.
*/ 

#include <util/global.h>
#include "CudaTestComposite.h"

#include <test/CompositeTestRunner.h>

using namespace Pscf;
using namespace Util;

int main(int argc, char* argv[])
{
   CudaTestComposite runner;

   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}

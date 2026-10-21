// Unit-test runner for the dvm-dos-tem daily-climate-mode changes.
// Build/run: `make unittest` from the repository root (see tests/README.md).
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "../../include/TEMLogger.h"
#include <iostream>

int main(int argc, char** argv) {
  // Model code logs through the global boost logger 'glg'; keep it quiet.
  std::streambuf* old = std::cout.rdbuf(nullptr);
  setup_logging("disabled", "all");
  std::cout.rdbuf(old);
  doctest::Context ctx;
  ctx.applyCommandLine(argc, argv);
  return ctx.run();
}

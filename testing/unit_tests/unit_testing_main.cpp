#include <gtest/gtest.h>
#include "logger/logger.hpp"
#include "testing/command_line.hpp"


int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

//  Logger::get().disable();

  // make command line arguments available to tests.
  std::vector<std::string> cmd_line(argv + 1, argv + argc);
  cmdline = cmd_line;

  return RUN_ALL_TESTS();
}

std::vector<std::string> cmdline;

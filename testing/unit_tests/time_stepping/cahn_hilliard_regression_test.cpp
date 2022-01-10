#include <gtest/gtest.h>

#include "initialization/puppeteer.hpp"
#include "logger/logger.hpp"
#include "testing/command_line.hpp"


TEST(CahnHilliardRegression, Step10)
{
  std::vector<std::string> cmd = cmdline;
  cmd.push_back("--max_time_steps=10");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  auto& res = Logger::get().getResidualMonitor();
  auto& sol = Logger::get().getSolutionMonitor();

  EXPECT_NEAR(7.32339e+01, res[0].second, 0.0001e1);   // res_c
  EXPECT_NEAR(-6.09467e-03, sol[0].second, 0.0001e-3); // min_c
  EXPECT_NEAR(5.89425e-03, sol[1].second, 0.0001e-3);  // max_c
  EXPECT_NEAR(2.46571e-01, sol[2].second, 0.0001e-1);  // grad_c_mag
}

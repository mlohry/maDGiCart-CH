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

  EXPECT_NEAR(73.233883471387173, res[0].second, 0.0001e1);   // res_c
  EXPECT_NEAR(-0.0060905932900116051, sol[0].second, 0.0001e-3); // min_c
  EXPECT_NEAR(0.0058983260761296171, sol[1].second, 0.0001e-3);  // max_c
  EXPECT_NEAR(0.24657124985861673, sol[2].second, 0.0001e-1);  // grad_c_mag
}

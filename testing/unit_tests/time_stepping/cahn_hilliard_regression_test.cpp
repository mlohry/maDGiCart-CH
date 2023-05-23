#include <gtest/gtest.h>

#include "initialization/puppeteer.hpp"
#include "logger/logger.hpp"
#include "testing/command_line.hpp"


TEST(CahnHilliardRegression, Step10RK3)
{
  std::vector<std::string> cmd = cmdline;
  cmd.push_back("--max_time_steps=10");
  cmd.push_back("--log_frequency=1");
  cmd.push_back("--time_step_size=1.e-6");
  cmd.push_back("--time_integrator=rk3ssp");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  auto& res = Logger::get().getResidualMonitor();
  auto& sol = Logger::get().getSolutionMonitor();

  EXPECT_NEAR(7.32339e+01, res[0].second, 0.0001e1);   // res_c
  EXPECT_NEAR(-6.09467e-03, sol[0].second, 0.0001e-3); // min_c
  EXPECT_NEAR(5.89425e-03, sol[1].second, 0.0001e-3);  // max_c
  EXPECT_NEAR(2.46571e-01, sol[2].second, 0.0001e-1);  // grad_c_mag
}


TEST(CahnHilliardRegression, Step10ODE23)
{
  std::vector<std::string> cmd = cmdline;
  cmd.push_back("--max_time_steps=10");
  cmd.push_back("--log_frequency=1");
  cmd.push_back("--time_step_size=1.e-6");
  cmd.push_back("--min_time_step_size=1.e-6");
  cmd.push_back("--max_time_step_size=1.e-6");
  cmd.push_back("--time_integrator=ode23");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  auto& res = Logger::get().getResidualMonitor();
  auto& sol = Logger::get().getSolutionMonitor();

  EXPECT_NEAR(7.32339e+01, res[0].second, 0.0001e1);   // res_c
  EXPECT_NEAR(-6.09467e-03, sol[0].second, 0.0001e-3); // min_c
  EXPECT_NEAR(5.89425e-03, sol[1].second, 0.0001e-3);  // max_c
  EXPECT_NEAR(2.46571e-01, sol[2].second, 0.0001e-1);  // grad_c_mag
}

TEST(CahnHilliardRegression, Step10PetscODE23)
{
  std::vector<std::string> cmd = cmdline;
  cmd.push_back("--max_time_steps=10");
  cmd.push_back("--log_frequency=1");
  cmd.push_back("--time_step_size=1.e-6");
  cmd.push_back("--min_time_step_size=1.e-6");
  cmd.push_back("--max_time_step_size=1.e-6");
  cmd.push_back("--time_integrator=petsc_ode23");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  auto& res = Logger::get().getResidualMonitor();
  auto& sol = Logger::get().getSolutionMonitor();

  EXPECT_NEAR(6.97340e+01, res[0].second, 0.0001e1);   // res_c
  EXPECT_NEAR(-6.09467e-03, sol[0].second, 0.0001e-3); // min_c
  EXPECT_NEAR(5.89425e-03, sol[1].second, 0.0001e-3);  // max_c
  EXPECT_NEAR(2.46571e-01, sol[2].second, 0.0001e-1);  // grad_c_mag
}

TEST(CahnHilliardRegression, DISABLED_PetscImplicitThin)
{
  const std::string config_file = std::string(SOURCEDIR) + "/configs/long_slow_thin_modified.dat";

  std::vector<std::string> cmd = cmdline;
  cmd.push_back("--max_time_steps=10");
  cmd.push_back("--config_file=" + config_file);
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  auto& res = Logger::get().getResidualMonitor();
  auto& sol = Logger::get().getSolutionMonitor();

  EXPECT_NEAR(638.90788002024851, res[0].second, 1.e-4);   // res_c
  EXPECT_NEAR(0.295477755185511, sol[0].second, 1.e-4); // min_c
  EXPECT_NEAR(0.30427039850611504, sol[1].second, 1.e-4);  // max_c
  EXPECT_NEAR(10.316461553738717, sol[2].second, 1.e-4);  // grad_c_mag
}

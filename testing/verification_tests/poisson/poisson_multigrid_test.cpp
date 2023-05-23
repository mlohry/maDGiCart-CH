#include <gtest/gtest.h>

#include "data_structures/exec_includes.hpp"
#include "initialization/puppeteer.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "testing/command_line.hpp"
#include "time_stepping/bdf_multigrid.hpp"


namespace {
double
poisson_l2error(const Puppeteer& puppeteer)
{
  const auto& time_integrator = puppeteer.timeIntegrator();
  const auto& geom            = dynamic_cast<const Discretization3DCart&>(puppeteer.spatialDiscretization());

  auto idx = read_access(geom.interiorIndices());
  auto sol = read_access(time_integrator.getCurrentSolutionState().getVec(0));
  auto x   = read_access(geom.x().asArray());
  auto y   = read_access(geom.y().asArray());
  auto z   = read_access(geom.z().asArray());

  ReduceSumReal squared_norm(0);

  maDGForAll(ii, 0, idx.size(), {
    const int     i        = idx[ii];
    const real_wp expected = sin(x[i]) * cos(y[i]) * cos(z[i]);
    squared_norm += pow(sol[i] - expected, 2.0);
  });
  return sqrt(squared_norm.get() / idx.size());
}
}  // namespace

TEST(PoissonMultigrid, CubeSingleGrid)
{
  std::vector<std::string> cmd         = cmdline;
  const std::string        config_file = std::string(SOURCEDIR) + "/configs/poisson_128cubed.dat";
  cmd.push_back("--config_file=" + config_file);
  cmd.push_back("--mg_levels=1");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  const auto& mg_solver = dynamic_cast<const BDFMultigrid&>(puppeteer.timeIntegrator());

  const auto ncycles = mg_solver.nCyclesToConverge();
  const auto rate    = mg_solver.convergenceRate();
  const auto l2err   = poisson_l2error(puppeteer);
  std::cout << "ncycles: " << ncycles << "\nrate: " << rate << "\nl2 error: " << l2err << "\n";

  EXPECT_EQ(ncycles, 2183);
  EXPECT_NEAR(rate, 0.993691, 0.001);
  EXPECT_NEAR(l2err, 4.8943269394794269e-05, 1.e-6);
}

TEST(PoissonMultigrid, Cube2V)
{
  std::vector<std::string> cmd         = cmdline;
  const std::string        config_file = std::string(SOURCEDIR) + "/configs/poisson_128cubed.dat";
  cmd.push_back("--config_file=" + config_file);
  cmd.push_back("--mg_levels=2");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  const auto& mg_solver = dynamic_cast<const BDFMultigrid&>(puppeteer.timeIntegrator());

  const auto ncycles = mg_solver.nCyclesToConverge();
  const auto rate    = mg_solver.convergenceRate();
  const auto l2err   = poisson_l2error(puppeteer);
  std::cout << "ncycles: " << ncycles << "\nrate: " << rate << "\nl2 error: " << l2err << "\n";

  EXPECT_EQ(ncycles, 477);
  EXPECT_NEAR(rate, 0.97142430398003676, 0.001);
  EXPECT_NEAR(l2err, 4.9087151442066223e-05, 1.e-6);
}

TEST(PoissonMultigrid, Cube3V)
{
  std::vector<std::string> cmd         = cmdline;
  const std::string        config_file = std::string(SOURCEDIR) + "/configs/poisson_128cubed.dat";
  cmd.push_back("--config_file=" + config_file);
  cmd.push_back("--mg_levels=3");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  const auto& mg_solver = dynamic_cast<const BDFMultigrid&>(puppeteer.timeIntegrator());

  const auto ncycles = mg_solver.nCyclesToConverge();
  const auto rate    = mg_solver.convergenceRate();
  const auto l2err   = poisson_l2error(puppeteer);
  std::cout << "ncycles: " << ncycles << "\nrate: " << rate << "\nl2 error: " << l2err << "\n";

  EXPECT_EQ(ncycles, 120);
  EXPECT_NEAR(rate, 0.89087274613233169, 0.001);
  EXPECT_NEAR(l2err, 5.0483117884540851e-05, 1.e-6);
}

TEST(PoissonMultigrid, Cube4V)
{
  std::vector<std::string> cmd         = cmdline;
  const std::string        config_file = std::string(SOURCEDIR) + "/configs/poisson_128cubed.dat";
  cmd.push_back("--config_file=" + config_file);
  cmd.push_back("--mg_levels=4");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  const auto& mg_solver = dynamic_cast<const BDFMultigrid&>(puppeteer.timeIntegrator());

  const auto ncycles = mg_solver.nCyclesToConverge();
  const auto rate    = mg_solver.convergenceRate();
  const auto l2err   = poisson_l2error(puppeteer);
  std::cout << "ncycles: " << ncycles << "\nrate: " << rate << "\nl2 error: " << l2err << "\n";

  EXPECT_EQ(ncycles, 38);
  EXPECT_NEAR(rate, 0.68947972738853158, 0.001);
  EXPECT_NEAR(l2err, 5.5462309342152126e-05, 1.e-6);
}

TEST(PoissonMultigrid, Cube5V)
{
  std::vector<std::string> cmd         = cmdline;
  const std::string        config_file = std::string(SOURCEDIR) + "/configs/poisson_128cubed.dat";
  cmd.push_back("--config_file=" + config_file);
  cmd.push_back("--mg_levels=5");
  cmd.push_back("--mg_nsmooth=1");
  cmd.push_back("--mg_interpolation=linear");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  const auto& mg_solver = dynamic_cast<const BDFMultigrid&>(puppeteer.timeIntegrator());

  const auto ncycles = mg_solver.nCyclesToConverge();
  const auto rate    = mg_solver.convergenceRate();
  const auto l2err   = poisson_l2error(puppeteer);
  std::cout << "ncycles: " << ncycles << "\nrate: " << rate << "\nl2 error: " << l2err << "\n";

  EXPECT_EQ(ncycles, 24);
  EXPECT_NEAR(rate, 0.55739672324878919, 0.1);
  EXPECT_NEAR(l2err, 5.62692e-05, 1.e-6);
}

TEST(PoissonMultigrid, Cube6V)
{
  std::vector<std::string> cmd         = cmdline;
  const std::string        config_file = std::string(SOURCEDIR) + "/configs/poisson_128cubed.dat";
  cmd.push_back("--config_file=" + config_file);
  cmd.push_back("--mg_levels=6");
  cmd.push_back("--mg_nsmooth=1");
  cmd.push_back("--mg_interpolation=linear");
  Puppeteer puppeteer(cmd);
  puppeteer.run();

  const auto& mg_solver = dynamic_cast<const BDFMultigrid&>(puppeteer.timeIntegrator());

  const auto ncycles = mg_solver.nCyclesToConverge();
  const auto rate    = mg_solver.convergenceRate();
  const auto l2err   = poisson_l2error(puppeteer);
  std::cout << "ncycles: " << ncycles << "\nrate: " << rate << "\nl2 error: " << l2err << "\n";

  EXPECT_EQ(ncycles, 20);
  EXPECT_NEAR(rate, 0.49943938499793233, 0.001);
  EXPECT_NEAR(l2err, 5.62692e-05, 5.e-6);
}

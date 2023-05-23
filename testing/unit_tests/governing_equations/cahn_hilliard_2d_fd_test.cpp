#include <gtest/gtest.h>

#include "typedefs.hpp"

#include "data_structures/exec_includes.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_2d_fd.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "testing/utils/order_of_accuracy.hpp"
#include "testing/utils/domain_generator.hpp"


TEST(CahnHilliard, RHSEvaluation)
{
  OrderOfAccuracy        accuracy_chrhs;
  const std::vector<int> mesh_sizes     = {16, 32, 64, 128, 256};
  const int              expected_order = 2;

  for (auto N : mesh_sizes) {
    Discretization2DCart   geom(generateTriplyPeriodicDomain(N));
    CahnHilliardParameters ch_params;
    CahnHilliard2DFD       ch_rhs(geom, ch_params);

    auto state     = ch_rhs.createSolutionState();
    auto dstate_dt = ch_rhs.createSolutionState();

    auto idx = read_access_host(geom.interiorIndices());
    auto f   = write_access_host(dynamic_cast<ScalarSolutionState2D&>(*state).c());
    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      f.getIJ(idx[ii], i, j);
      f(i, j) = sin(x(i, j)) * cos(y(i, j));
    });

    ch_rhs.evalRHSImpl(*state, 0, *dstate_dt);

    ReduceSumRealHostSeq squared_norm(0);
    const real_wp m    = ch_params.m();
    const real_wp eps2 = ch_params.eps2();
    const real_wp sigma = ch_params.sigma();

    auto rhs = read_access_host(dynamic_cast<ScalarSolutionState2D&>(*dstate_dt).c());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      f.getIJ(idx[ii], i, j);

      const real_wp sinx = sin(x(i, j));
      const real_wp siny = sin(y(i, j));
      const real_wp cosx = cos(x(i, j));
      const real_wp cosy = cos(y(i, j));

      // laplacian (c^3 -c)
      const real_wp expected_laplacian =
          2 * sinx * cosy *
          (3 * sinx * sinx * siny * siny + 3 * cosx * cosx * cosy * cosy - 3 * sinx * sinx * cosy * cosy + 1);

      // biharmonic(c)
      const real_wp expected_biharmonic = 4 * sinx * cosy;

      const real_wp expected_linear = sinx * cosy - m;

      const real_wp expected = -eps2*expected_biharmonic + expected_laplacian - sigma * expected_linear;

      squared_norm += pow(expected - rhs(i, j), 2.0);
    });

    const real_wp l2error = std::sqrt(squared_norm.get() / geom.nInteriorPoints());

    accuracy_chrhs.addSolutionError(expected_order, N, l2error);
  }

  std::cout << "cartesian finite difference Cahn-Hilliard error:\n" << accuracy_chrhs << "\n";
  EXPECT_NEAR(accuracy_chrhs.getConvergenceRate().at(expected_order).back(), expected_order, 0.01);
}

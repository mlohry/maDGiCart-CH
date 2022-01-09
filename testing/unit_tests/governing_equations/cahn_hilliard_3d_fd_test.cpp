#include <gtest/gtest.h>

#include "typedefs.hpp"

#include "data_structures/exec_includes.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_3d_fd.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "testing/utils/order_of_accuracy.hpp"

#include "file_io/vtk_solution_writer.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_initial_conditions.hpp"
#include "governing_equations/residual_calculator.hpp"
#include "program_options/program_options_parser.hpp"
#include "time_stepping/rk3ssp.hpp"
#include "testing/command_line.hpp"

TEST(CahnHilliard3D, RHSEvaluation)
{
  OrderOfAccuracy        accuracy_chrhs;
  const std::vector<int> mesh_sizes     = {16, 32, 64, 128};
  const int              expected_order = 2;

  for (auto N : mesh_sizes) {
    Discretization3DCart   geom(N, 2, -M_PI, M_PI, -M_PI, -M_PI);
    CahnHilliardParameters ch_params;
    CahnHilliard3DFD       ch_rhs(geom, ch_params);

    auto state     = ch_rhs.createSolutionState();
    auto dstate_dt = ch_rhs.createSolutionState();

    auto idx = read_access_host(geom.interiorIndices());
    auto f   = write_access_host(dynamic_cast<CahnHilliardState3D&>(*state).c());
    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());
    auto z   = read_access_host(geom.z());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      f.getIJK(idx[ii], i, j, k);
      f(i, j, k) = sin(x(i, j, k)) * cos(y(i, j, k)) * cos(z(i, j, k));
    });

    ch_rhs.evalRHSImpl(*state, 0, *dstate_dt);

    ReduceSumRealHostSeq squared_norm(0);
    const real_wp        m     = ch_params.m();
    const real_wp        eps2  = ch_params.eps2();
    const real_wp        sigma = ch_params.sigma();

    auto rhs = read_access_host(dynamic_cast<CahnHilliardState3D&>(*dstate_dt).c());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      f.getIJK(idx[ii], i, j, k);
      const real_wp sinx = sin(x(i, j, k));
      const real_wp siny = sin(y(i, j, k));
      const real_wp sinz = sin(z(i, j, k));
      const real_wp cosx = cos(x(i, j, k));
      const real_wp cosy = cos(y(i, j, k));
      const real_wp cosz = cos(z(i, j, k));

      // laplacian (c^3 -c)
      // clang-format off
      const real_wp expected_laplacian_c3 =
          - real_wp(9) * pow(sinx,3.0) * pow(cosy,3.0) * pow(cosz,3.0)
          + real_wp(6) * pow(sinx,3.0) * pow(siny,2.0) * cosy * pow(cosz,3.0)
          + real_wp(6) * pow(sinx,3.0) * pow(cosy, 3.0) * pow(sinz,2.0) * cosz
          + real_wp(6) * sinx * pow(cosx,2.0) * pow(cosy,3.0) * pow(cosz,3.0);

      const real_wp expected_laplacian_c =
          - real_wp(3) * sinx * cosy * cosz;

      const real_wp expected_biharmonic = 9.0 * sinx * cosy * cosz;

      const real_wp expected_linear = sinx*cosy*cosz - m;

      const real_wp expected = -eps2 * expected_biharmonic + (expected_laplacian_c3-expected_laplacian_c) - sigma *  expected_linear;

      squared_norm += pow(expected - rhs(i, j, k), 2.0);

      // clang-format on
    });

    const real_wp l2error = std::sqrt(squared_norm.get() / geom.nInteriorPoints());

    accuracy_chrhs.addSolutionError(expected_order, N, l2error);
  }

  std::cout << "cartesian finite difference Cahn-Hilliard error:\n" << accuracy_chrhs << "\n";
  EXPECT_NEAR(accuracy_chrhs.getConvergenceRate().at(expected_order).back(), expected_order, 0.01);
}

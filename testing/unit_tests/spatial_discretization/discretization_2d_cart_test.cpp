#include <gtest/gtest.h>

#include "typedefs.hpp"

#include "data_structures/exec_includes.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "testing/utils/order_of_accuracy.hpp"


TEST(Discretization2DCart, Laplacian)
{
  OrderOfAccuracy        accuracy_laplacian;
  const std::vector<int> mesh_sizes     = {16, 32, 64, 128, 256};
  const int              expected_order = 2;

  for (auto N : mesh_sizes) {
    Discretization2DCart geom(N, 2, -M_PI, M_PI, -M_PI);
    auto state = geom.createRealArray();
    auto del2  = geom.createRealArray();

    auto idx = read_access_host(geom.interiorIndices());
    auto f   = write_access_host((*state));
    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      f.getIJ(idx[ii], i, j);
      f(i, j) = sin(x(i, j)) * cos(y(i, j));
    });

    geom.applyPeriodicBoundaryConditions(*state);
    geom.laplacian(*state, *del2);

    ReduceSumRealHostSeq squared_norm(0);
    auto          laplacian = read_access_host((*del2));

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      f.getIJ(idx[ii], i, j);
      const real_wp expected = -2.0 * sin(x(i, j)) * cos(y(i, j));
      squared_norm += pow(expected - laplacian(i, j), 2.0);
    });

    const real_wp l2error = std::sqrt(squared_norm.get() / geom.nInteriorPoints());

    accuracy_laplacian.addSolutionError(expected_order, N, l2error);
  }

  std::cout << "cartesian finite difference laplacian error:\n" << accuracy_laplacian << "\n";
  EXPECT_NEAR(accuracy_laplacian.getConvergenceRate().at(expected_order).back(), expected_order, 0.01);
}


TEST(Discretization2DCart, Biharmonic)
{
  OrderOfAccuracy        accuracy_biharmonic;
  const std::vector<int> mesh_sizes     = {16, 32, 64, 128, 256};
  const int              expected_order = 2;

  for (auto N : mesh_sizes) {
    Discretization2DCart geom(N, 2, -M_PI, M_PI, -M_PI);
    auto state = geom.createRealArray();
    auto del4  = geom.createRealArray();

    auto idx = read_access_host(geom.interiorIndices());
    auto f   = write_access_host((*state));
    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      f.getIJ(idx[ii], i, j);
      f(i, j) = sin(x(i, j)) * cos(y(i, j));
    });
    geom.applyPeriodicBoundaryConditions(*state);
    geom.biharmonic(*state, *del4);

    ReduceSumRealHostSeq squared_norm(0);

    auto biharmonic = read_access_host((*del4));

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      f.getIJ(idx[ii], i, j);

      const real_wp expected = 4.0 * sin(x(i, j)) * cos(y(i, j));

      squared_norm += pow(expected - biharmonic(i, j), 2.0);
    });

    const real_wp l2error = std::sqrt(squared_norm.get() / geom.nInteriorPoints());

    accuracy_biharmonic.addSolutionError(expected_order, N, l2error);
  }

  std::cout << "cartesian finite difference biharmonic error:\n" << accuracy_biharmonic << "\n";
  EXPECT_NEAR(accuracy_biharmonic.getConvergenceRate().at(expected_order).back(), expected_order, 0.01);
}

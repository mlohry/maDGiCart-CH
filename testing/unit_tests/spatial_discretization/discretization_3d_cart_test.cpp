#include <gtest/gtest.h>

#include "typedefs.hpp"

#include "data_structures/exec_includes.hpp"
#include "data_structures/managed_array_3d.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "testing/utils/order_of_accuracy.hpp"


TEST(Discretiztion3DCart, DataStructures)
{
  const int                    N     = 2;
  const int                    Nhalo = 1;
  ManagedArray3DOwning<double> array3d({}, "arrayname", N, N, N, Nhalo);

  auto arr = read_access_host(array3d);

  int idx1d = 0;
  for (int i = -Nhalo; i < N + Nhalo; ++i) {
    for (int j = -Nhalo; j < N + Nhalo; ++j) {
      for (int k = -Nhalo; k < N + Nhalo; ++k) {
        int ii, jj, kk;
        arr.getIJK(idx1d, ii, jj, kk);
        ASSERT_EQ(i, ii);
        ASSERT_EQ(j, jj);
        ASSERT_EQ(k, kk);
        idx1d++;
      }
    }
  }
}


TEST(Discretiztion3DCart, Laplacian)
{
  OrderOfAccuracy        accuracy_laplacian;
  const std::vector<int> mesh_sizes     = {16, 32, 64, 128};
  const int              expected_order = 2;

  for (auto N : mesh_sizes) {
    Discretization3DCart geom(N, 2, -M_PI, M_PI, -M_PI, -M_PI);
    auto                 state = geom.createRealArray();
    auto                 del2  = geom.createRealArray();

    auto idx = read_access_host(geom.interiorIndices());
    auto f   = write_access_host((*state));
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

    geom.applyPeriodicBoundaryConditions(*state);
    geom.laplacian(*state, *del2);

    ReduceSumRealHostSeq squared_norm(0);
    auto                 laplacian = read_access_host((*del2));

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      f.getIJK(idx[ii], i, j, k);
      const real_wp expected = -3.0 * sin(x(i, j, k)) * cos(y(i, j, k)) * cos(z(i, j, k));
      squared_norm += pow(expected - laplacian(i, j, k), 2.0);
    });

    const real_wp l2error = std::sqrt(squared_norm.get() / geom.nInteriorPoints());
    accuracy_laplacian.addSolutionError(expected_order, N, l2error);
  }

  std::cout << "cartesian finite difference laplacian error:\n" << accuracy_laplacian << "\n";
  EXPECT_NEAR(accuracy_laplacian.getConvergenceRate().at(expected_order).back(), expected_order, 0.01);
}

TEST(Discretiztion3DCart, Biharmonic)
{
  OrderOfAccuracy        accuracy_biharmonic;
  const std::vector<int> mesh_sizes     = {16, 32, 64, 128};
  const int              expected_order = 2;

  for (auto N : mesh_sizes) {
    Discretization3DCart geom(N, 2, -M_PI, M_PI, -M_PI, -M_PI);
    auto                 state = geom.createRealArray();
    auto                 del4  = geom.createRealArray();

    auto idx = read_access_host(geom.interiorIndices());
    auto f   = write_access_host((*state));
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

    geom.applyPeriodicBoundaryConditions(*state);
    geom.biharmonic(*state, *del4);

    ReduceSumRealHostSeq squared_norm(0);

    auto biharmonic = read_access_host((*del4));

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      f.getIJK(idx[ii], i, j, k);

      const real_wp expected = 9.0 * sin(x(i, j, k)) * cos(y(i, j, k)) * cos(z(i, j, k));

      squared_norm += pow(expected - biharmonic(i, j, k), 2.0);
//      std::cout << "i,j,k: " << i << " " << j << " "<< k << " " << biharmonic(i,j,k) << " expected " << expected << "\n";
    });

    const real_wp l2error = std::sqrt(squared_norm.get() / geom.nInteriorPoints());

    accuracy_biharmonic.addSolutionError(expected_order, N, l2error);
  }
  std::cout << "cartesian finite difference biharmonic error:\n" << accuracy_biharmonic << "\n";
  EXPECT_NEAR(accuracy_biharmonic.getConvergenceRate().at(expected_order).back(), expected_order, 0.01);
}
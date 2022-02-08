#include <gtest/gtest.h>

#include "typedefs.hpp"

#include "data_structures/exec_includes.hpp"
#include "file_io/vtk_solution_reader.hpp"
#include "file_io/vtk_solution_writer.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_2d_fd.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_3d_fd.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"


TEST(FileIO, Discretization2DCartWriteAndRead)
{
  const std::string filename_without_ext = "fileio_unittest";
  const std::string filename_with_ext    = filename_without_ext + ".vts";

  const int    N     = 16;
  const int    Nhalo = 2;
  const double xbeg  = -2;
  const double xend  = 3;
  const double ybeg  = -1;

  CartesianDomainDefinition domain;
  domain.nx = domain.ny = domain.nz = N;

  domain.xbeg  = xbeg;
  domain.xend  = xend;
  domain.ybeg  = ybeg;
  domain.nhalo = Nhalo;


  Discretization2DCart   geom(domain);
  CahnHilliardParameters ch_params;
  CahnHilliard2DFD       rhs(geom, ch_params);

  //  write a specific solution file to disk
  {
    auto state = rhs.createSolutionState();


    auto c = write_access_host(dynamic_cast<CahnHilliardState&>(*state).c());

    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());
    auto idx = read_access_host(geom.interiorIndices());


    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      c.getIJ(idx[ii], i, j);
      c(i, j) = sin(x(i, j)) * cos(y(i, j));
    });

    write_solution_to_vtk(*state, geom, filename_without_ext);
  }


  // read back the solution file
  {
    VTKSolutionReader file_reader(filename_with_ext);

    const auto domain_defn = file_reader.getCartesianDomain();
    EXPECT_EQ(domain_defn.nx, N);
    EXPECT_EQ(domain_defn.ny, N);
    EXPECT_EQ(domain_defn.nz, 0);
    EXPECT_EQ(domain_defn.xbeg, xbeg);
    EXPECT_EQ(domain_defn.xend, xend);
    EXPECT_EQ(domain_defn.ybeg, ybeg);
    EXPECT_EQ(domain_defn.zbeg, 0);

    auto input_state = rhs.createSolutionState();
    file_reader.setSolution(geom.interiorIndices(), *input_state);


    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());
    auto idx = read_access_host(geom.interiorIndices());
    auto c   = read_access_host(dynamic_cast<CahnHilliardState&>(*input_state).c());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      c.getIJ(idx[ii], i, j);
      EXPECT_NEAR(c(i, j), sin(x(i, j)) * cos(y(i, j)), 1.e-16);
    });
  }
}


TEST(FileIO, Discretization3DCartWriteAndRead)
{
  const std::string filename_without_ext = "fileio3d_unittest";
  const std::string filename_with_ext    = filename_without_ext + ".vts";

  const int    N     = 16;
  const int    Nhalo = 2;
  const double xbeg  = -2;
  const double xend  = 3;
  const double ybeg  = -1;
  const double zbeg  = 1;

  CartesianDomainDefinition domain;
  domain.nx = domain.ny = domain.nz = N;

  domain.xbeg  = xbeg;
  domain.xend  = xend;
  domain.ybeg  = ybeg;
  domain.zbeg  = zbeg;
  domain.nhalo = Nhalo;

  Discretization3DCart   geom(domain);
  CahnHilliardParameters ch_params;
  CahnHilliard3DFD       rhs(geom, ch_params);

  //  write a specific solution file to disk
  {
    auto state = rhs.createSolutionState();


    auto c = write_access_host(dynamic_cast<CahnHilliardState3D&>(*state).c());

    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());
    auto z   = read_access_host(geom.z());
    auto idx = read_access_host(geom.interiorIndices());


    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      c.getIJK(idx[ii], i, j, k);
      c(i, j, k) = sin(x(i, j, k)) * cos(y(i, j, k)) * cos(z(i, j, k));
    });

    write_solution_to_vtk(*state, geom, filename_without_ext);
  }


  // read back the solution file
  {
    VTKSolutionReader file_reader(filename_with_ext);

    const auto domain_defn = file_reader.getCartesianDomain();
    EXPECT_EQ(domain_defn.nx, N);
    EXPECT_EQ(domain_defn.ny, N);
    EXPECT_EQ(domain_defn.nz, N);
    EXPECT_EQ(domain_defn.xbeg, xbeg);
    EXPECT_EQ(domain_defn.xend, xend);
    EXPECT_EQ(domain_defn.ybeg, ybeg);
    EXPECT_EQ(domain_defn.zbeg, zbeg);

    auto input_state = rhs.createSolutionState();
    file_reader.setSolution(geom.interiorIndices(), *input_state);


    auto x   = read_access_host(geom.x());
    auto y   = read_access_host(geom.y());
    auto z   = read_access_host(geom.z());
    auto idx = read_access_host(geom.interiorIndices());
    auto c   = read_access_host(dynamic_cast<CahnHilliardState3D&>(*input_state).c());

    maDGForAllHost(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      c.getIJK(idx[ii], i, j, k);
      EXPECT_NEAR(c(i, j, k), sin(x(i, j, k)) * cos(y(i, j, k)) * cos(z(i, j, k)), 1.e-16);
    });
  }
}

#include "vtk_solution_writer.hpp"

#include <fstream>
#include <iomanip>
#include <limits>

#include "data_structures/solution_state.hpp"
#include "logger/logger.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_state.hpp"


void
write_solution_to_vtk(const SolutionState& state, const Discretization2DCart& geom, std::string filename_withoutext)
{
  const std::string filename = filename_withoutext + ".vts";
  auto              timer    = Logger::get().timer("Writing VTK file " + filename);
  std::ofstream     fout;
  fout.open(filename);


  const CahnHilliardState& cstate = dynamic_cast<const CahnHilliardState&>(state);

  const int ni  = geom.ni();
  auto      x   = read_access_host(geom.x());
  auto      y   = read_access_host(geom.y());
  auto      val = read_access_host(cstate.c());

  fout << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
  fout << "  <StructuredGrid WholeExtent=\"0 " << ni - 1 << " 0 " << ni - 1 << " 0 0\">\n";
  fout << "    <Piece Extent=\" 0 " << ni - 1 << " 0 " << ni - 1 << " 0 0\">\n";
  fout << "      <PointData>\n";
  fout << "      </PointData>\n";
  fout << "      <CellData>\n";

  fout << "    <DataArray type=\"Float64\" Name=\"c\" format=\"ascii\">\n";
  for (int i = 0; i < ni - 1; ++i) {
    for (int j = 0; j < ni - 1; ++j) {
      fout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << val(i, j) << "\n";
    }
  }
  fout << "\n    </DataArray>\n";
  fout << "      </CellData>\n";
  fout << "      <Points>\n";
  fout << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" "
          "format=\"ascii\">\n";

  for (int i = 0; i < ni; ++i) {
    for (int j = 0; j < ni; ++j) {
      fout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << x(i, j) << " " << y(i, j) << " 0.0\n";
    }
  }

  fout << "        </DataArray>\n";
  fout << "      </Points>\n";
  fout << "    </Piece>\n";
  fout << "  </StructuredGrid>\n";
  fout << "</VTKFile>\n";
}

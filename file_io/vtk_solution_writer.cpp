#include "vtk_solution_writer.hpp"

#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <limits>

#include "data_structures/scalar_solution_state.hpp"
#include "data_structures/solution_state.hpp"
#include "logger/logger.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "time_stepping/time_integrator.hpp"


void
write_solution_to_vtk(
    const std::vector<std::pair<std::string, double>>& iteration_status,
    const SolutionState&                               state,
    const SpatialDiscretization&                       geom,
    std::string                                        filename_withoutext)
{
  try {
    const auto& g2d = dynamic_cast<const Discretization2DCart&>(geom);
    write_solution_to_vtk_2d(iteration_status, state, g2d, filename_withoutext);
  }
  catch (...) {
    const auto& g3d = dynamic_cast<const Discretization3DCart&>(geom);
    write_solution_to_vtk_3d(iteration_status, state, g3d, filename_withoutext);
  }
}


void
write_solution_to_vtk_2d(
    const std::vector<std::pair<std::string, double>>& iteration_status,
    const SolutionState&                               state,
    const Discretization2DCart&                        geom,
    std::string                                        filename_withoutext)
{
  const int                     iterwidth = int(log10(Options::get().max_time_steps()) + 1);
  std::map<std::string, double> iter_map;
  std::copy(iteration_status.begin(), iteration_status.end(), std::inserter(iter_map, iter_map.begin()));

  std::ostringstream ss;
  ss << filename_withoutext;
  ss << "_step_" << std::setw(iterwidth) << std::setfill('0') << iter_map.at("iter");
  ss << ".vts";

  const std::string filename = ss.str();

  auto          timer = Logger::get().timer("Writing VTK file " + filename);
  std::ofstream fout;
  fout.open(filename);
  fout.imbue(std::locale("C")); // prevent commas in integer output

  const ScalarSolutionState2D& cstate = dynamic_cast<const ScalarSolutionState2D&>(state);

  const int ni = geom.ni();

  auto x   = read_access_host(geom.xvertex());
  auto y   = read_access_host(geom.yvertex());
  auto val = read_access_host(cstate.c());

  fout << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
  fout << "  <StructuredGrid WholeExtent=\"0 " << ni << " 0 " << ni << " 0 0\">\n";
  fout << "    <FieldData>\n";
  fout << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">\n";
  fout << std::scientific << iter_map.at("time") << "\n";
  fout << "     </DataArray>\n";
  fout << "     <DataArray type=\"Int32\" Name=\"IterationValue\" NumberOfTuples=\"1\">\n";
  fout << int(iter_map.at("iter")) << "\n";
  fout << "     </DataArray>\n";
  fout << "     <DataArray type=\"Float64\" Name=\"TimeStepSize\" NumberOfTuples=\"1\">\n";
  fout << std::scientific << iter_map.at("dt") << "\n";
  fout << "     </DataArray>\n";
  fout << "    </FieldData>\n";
  fout << "    <Piece Extent=\" 0 " << ni << " 0 " << ni << " 0 0\">\n";
  fout << "      <PointData>\n";
  fout << "      </PointData>\n";
  fout << "      <CellData>\n";

  fout << "    <DataArray type=\"Float64\" Name=\"c\" format=\"ascii\">\n";
  for (int j = 0; j < ni; ++j) {
    for (int i = 0; i < ni; ++i) {
      fout << std::setprecision(std::numeric_limits<real_wp>::digits10 + 1) << val(i, j) << "\n";
    }
  }
  fout << "\n    </DataArray>\n";
  fout << "      </CellData>\n";
  fout << "      <Points>\n";
  fout << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" "
          "format=\"ascii\">\n";

  for (int j = 0; j < ni + 1; ++j) {
    for (int i = 0; i < ni + 1; ++i) {
      fout << std::setprecision(std::numeric_limits<real_wp>::digits10 + 1) << x(i, j) << " " << y(i, j) << " 0.0\n";
    }
  }

  fout << "        </DataArray>\n";
  fout << "      </Points>\n";
  fout << "    </Piece>\n";
  fout << "  </StructuredGrid>\n";
  fout << "</VTKFile>\n";
}


void
write_solution_to_vtk_3d(
    const std::vector<std::pair<std::string, double>>& iteration_status,
    const SolutionState&                               state,
    const Discretization3DCart&                        geom,
    std::string                                        filename_withoutext)
{
  const int                     iterwidth = std::max(int(log10(Options::get().max_time_steps()) + 1), 4);
  std::map<std::string, double> iter_map;
  std::copy(iteration_status.begin(), iteration_status.end(), std::inserter(iter_map, iter_map.begin()));

  std::string iterstring = boost::lexical_cast<std::string>(iter_map.at("iter"));
  iterstring.erase(std::remove(iterstring.begin(), iterstring.end(), ','), iterstring.end());


  std::ostringstream ss;
  ss << filename_withoutext;
  ss << "_step_" << std::setw(iterwidth) << std::setfill('0') << iterstring;
  ss << ".vts";

  const std::string filename = ss.str();
  auto              timer    = Logger::get().timer("Writing VTK file " + filename);
  std::ofstream     fout;
  fout.open(filename);
  fout.imbue(std::locale("C")); // prevent commas in integer output

  const ScalarSolutionState3D& cstate = dynamic_cast<const ScalarSolutionState3D&>(state);

  const int ni  = geom.ni();
  const int nj  = geom.nj();
  const int nk  = geom.nk();
  auto      x   = read_access_host(geom.xvertex());
  auto      y   = read_access_host(geom.yvertex());
  auto      z   = read_access_host(geom.zvertex());
  auto      val = read_access_host(cstate.c());


  fout << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
  fout << "  <StructuredGrid WholeExtent=\"0 " << ni << " 0 " << nj << " 0 " << nk << "\">\n";
  fout << "    <FieldData>\n";
  fout << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">\n";
  fout << std::scientific << iter_map.at("time") << "\n";
  fout << "     </DataArray>\n";
  fout << "     <DataArray type=\"Int32\" Name=\"IterationValue\" NumberOfTuples=\"1\">\n";
  fout << int(iter_map.at("iter")) << "\n";
  fout << "     </DataArray>\n";
  fout << "     <DataArray type=\"Float64\" Name=\"TimeStepSize\" NumberOfTuples=\"1\">\n";
  fout << std::scientific << iter_map.at("dt") << "\n";
  fout << "     </DataArray>\n";
  fout << "    </FieldData>\n";
  fout << "    <Piece Extent=\" 0 " << ni << " 0 " << nj << " 0 " << nk << "\">\n";
  fout << "      <PointData>\n";
  fout << "      </PointData>\n";
  fout << "      <CellData>\n";
  fout << "    <DataArray type=\"Float64\" Name=\"c\" format=\"ascii\">\n";
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        fout << std::setprecision(std::numeric_limits<real_wp>::digits10 + 1) << val(i, j, k) << "\n";
      }
    }
  }
  fout << "\n    </DataArray>\n";
  fout << "      </CellData>\n";
  fout << "      <Points>\n";
  fout << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" "
          "format=\"ascii\">\n";

  for (int k = 0; k < nk + 1; ++k) {
    for (int j = 0; j < nj + 1; ++j) {
      for (int i = 0; i < ni + 1; ++i) {
        fout << std::setprecision(std::numeric_limits<real_wp>::digits10 + 1) << x(i, j, k) << " " << y(i, j, k) << " "
             << z(i, j, k) << "\n";
      }
    }
  }

  fout << "        </DataArray>\n";
  fout << "      </Points>\n";
  fout << "    </Piece>\n";
  fout << "  </StructuredGrid>\n";
  fout << "</VTKFile>\n";
}

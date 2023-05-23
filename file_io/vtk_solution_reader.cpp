#include "vtk_solution_reader.hpp"

#include "data_structures/scalar_solution_state.hpp"
#include "data_structures/solution_state.hpp"
#include "logger/logger.hpp"
#include "time_stepping/time_integrator_options.hpp"
#include "tinyxml2.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>


CartesianDomainDefinition
VTKSolutionReader::getCartesianDomain()
{
  tinyxml2::XMLDocument doc;
  auto                  timer = Logger::get().timer("VTKSolutionReader::getCartesianDomain() reading " + filename_);

  doc.LoadFile(filename_.c_str());

  const char* value_attribute = doc.FirstChildElement("VTKFile")
                                    ->FirstChildElement("StructuredGrid")
                                    ->FirstChildElement("Piece")
                                    ->FirstAttribute()
                                    ->Value();
  std::string attribute_string(value_attribute);
  // remove unwanted leading or trailing whitespace
  boost::trim_if(attribute_string, boost::is_any_of("\t "));


  // split by spsaces
  std::vector<std::string> tokens;
  boost::algorithm::split(tokens, attribute_string, boost::is_any_of("\t "), boost::token_compress_on);

  // cast to ints, will throw on failure
  const int nxbeg = boost::lexical_cast<int>(tokens[0]);
  const int nxend = boost::lexical_cast<int>(tokens[1]);
  const int nybeg = boost::lexical_cast<int>(tokens[2]);
  const int nyend = boost::lexical_cast<int>(tokens[3]);
  const int nzbeg = boost::lexical_cast<int>(tokens[4]);
  const int nzend = boost::lexical_cast<int>(tokens[5]);

  Logger::get().FatalAssert(nxbeg == 0, "VTKSolutionReader nxbeg == 0");
  Logger::get().FatalAssert(nybeg == 0, "VTKSolutionReader nybeg == 0");
  Logger::get().FatalAssert(nzbeg == 0, "VTKSolutionReader nzbeg == 0");

  CartesianDomainDefinition domain;
  domain.nx  = nxend;
  domain.ny  = nyend;
  dimension_ = 2;
  if (nzend == nzbeg) {
    domain.nz = 0;
  }
  else {
    domain.nz  = nzend;
    dimension_ = 3;
  }

  auto dataarray = doc.FirstChildElement("VTKFile")
                       ->FirstChildElement("StructuredGrid")
                       ->FirstChildElement("Piece")
                       ->FirstChildElement("CellData")
                       ->FirstChildElement("DataArray")
                       ->GetText();

  //  printf("Points GetText: %s\n", dataarray);
  std::string cell_data_string(dataarray);
  // remove unwanted leading or trailing whitespace
  boost::trim_if(cell_data_string, boost::is_any_of("\t \n"));

  std::vector<std::string> cell_data_tokens;
  boost::algorithm::split(cell_data_tokens, cell_data_string, boost::is_any_of("\t\n "), boost::token_compress_on);

  //  for (auto it = c
  auto points_array = doc.FirstChildElement("VTKFile")
                          ->FirstChildElement("StructuredGrid")
                          ->FirstChildElement("Piece")
                          ->FirstChildElement("Points")
                          ->FirstChildElement("DataArray")
                          ->GetText();

  std::string points_data_string(points_array);
  boost::trim_if(points_data_string, boost::is_any_of("\t \n"));
  std::vector<std::string> point_data_tokens;
  boost::algorithm::split(point_data_tokens, points_data_string, boost::is_any_of("\t\n "), boost::token_compress_on);


  const int npoints = point_data_tokens.size() / 3;
  if (dimension_ == 2) {
    Logger::get().FatalAssert(
        npoints == (nxend + 1) * (nyend + 1),
        "VTKSolutionReader npoints == (nx+1)*(ny+1), got npoints: " + std::to_string(npoints) +
            " nx: " + std::to_string(nxend) + " ny: " + std::to_string(nyend));
  }
  if (dimension_ == 3) {
    Logger::get().FatalAssert(
        npoints == (nxend + 1) * (nyend + 1) * (nzend + 1), "VTKSolutionReader npoints == (nx+1)*(ny+1)*(nz+1)");
  }


  std::vector<double> x, y, z;
  for (size_t i = 0; i < point_data_tokens.size(); i += 3) {
    x.push_back(boost::lexical_cast<double>(point_data_tokens[i]));
    y.push_back(boost::lexical_cast<double>(point_data_tokens[i + 1]));
    z.push_back(boost::lexical_cast<double>(point_data_tokens[i + 2]));
  }

  domain.xbeg = *min_element(x.begin(), x.end());
  domain.xend = *max_element(x.begin(), x.end());
  domain.ybeg = *min_element(y.begin(), y.end());
  domain.zbeg = *min_element(z.begin(), z.end());

  return domain;
}


void
VTKSolutionReader::setSolution(const IndexArray& interior_indices, SolutionState& state)
{
  tinyxml2::XMLDocument doc;
  auto                  timer = Logger::get().timer("VTKSolutionReader::setSolution() reading " + filename_);
  doc.LoadFile(filename_.c_str());

  auto solnarray = doc.FirstChildElement("VTKFile")
                       ->FirstChildElement("StructuredGrid")
                       ->FirstChildElement("Piece")
                       ->FirstChildElement("CellData")
                       ->FirstChildElement("DataArray")
                       ->GetText();

  std::string cell_data_string(solnarray);
  // remove unwanted leading or trailing whitespace
  boost::trim_if(cell_data_string, boost::is_any_of("\t \n"));

  std::vector<std::string> cell_data_tokens;
  boost::algorithm::split(cell_data_tokens, cell_data_string, boost::is_any_of("\t\n "), boost::token_compress_on);


  std::vector<double> solution(cell_data_tokens.size());
  for (size_t i = 0; i < solution.size(); ++i) {
    solution[i] = boost::lexical_cast<double>(cell_data_tokens[i]);
  }

  const auto domain = getCartesianDomain();

  if (!domain.nz) {  // 2d
    ScalarSolutionState2D& state2d = dynamic_cast<ScalarSolutionState2D&>(state);
    auto               f0      = write_access_host(state2d.c());
    int                solidx  = 0;
    for (int j = 0; j < domain.ny; ++j) {
      for (int i = 0; i < domain.nx; ++i) {
        f0(i, j) = solution[solidx];
        solidx++;
      }
    }
  }
  else {
    ScalarSolutionState3D& state3d = dynamic_cast<ScalarSolutionState3D&>(state);
    auto                 f0      = write_access_host(state3d.c());
    int                  solidx  = 0;
    for (int k = 0; k < domain.nz; ++k) {
      for (int j = 0; j < domain.ny; ++j) {
        for (int i = 0; i < domain.nx; ++i) {
          f0(i, j, k) = solution[solidx];
          solidx++;
        }
      }
    }
  }
}


void
VTKSolutionReader::setInitialTimeStep(TimeIntegratorOptions& opts)
{
  tinyxml2::XMLDocument doc;
  auto                  timer = Logger::get().timer("VTKSolutionReader::setSolution() reading " + filename_);
  doc.LoadFile(filename_.c_str());

  auto field_data =
      doc.FirstChildElement("VTKFile")->FirstChildElement("StructuredGrid")->FirstChildElement("FieldData");
  if (field_data) {
    for (auto* e = field_data->FirstChildElement("DataArray"); e != NULL; e = e->NextSiblingElement("DataArray")) {
      auto name_str = e->Attribute("Name");
      if (name_str) {
        auto name = std::string(name_str);
        Logger::get().InfoMessage("Read DataArray name " + name);
        auto        text_data = e->GetText();
        std::string data_string(text_data);
        // remove unwanted leading or trailing whitespace
        boost::trim_if(data_string, boost::is_any_of("\t \n"));
        std::vector<std::string> cell_data_tokens;
        boost::algorithm::split(cell_data_tokens, data_string, boost::is_any_of("\t\n "), boost::token_compress_on);

        if (name == "TimeValue") {
          opts.t0_ = boost::lexical_cast<real_t>(cell_data_tokens[0]);
          Logger::get().InfoMessage("Restarting at time " + std::to_string(opts.t0_));
        }
        else if (name == "IterationValue") {
          opts.initial_step_ = boost::lexical_cast<int>(cell_data_tokens[0]);
          Logger::get().InfoMessage("Restarting at iteration number " + std::to_string(opts.initial_step_));
        }
        else if (name == "TimeStepSize") {
          opts.dt_initial_ = boost::lexical_cast<real_t>(cell_data_tokens[0]);
          Logger::get().InfoMessage("Restarting with time step size " + std::to_string(opts.dt_initial_));
        }
      }
    }
  }
}
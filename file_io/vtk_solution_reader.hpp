#pragma once


#include "spatial_discretization/spatial_discretization.hpp"
#include "data_structures/managed_array.hpp"


class SolutionState;


class VTKSolutionReader /* : inherit public InitialConditions? better, move that into /initialization */ {
 public:
  VTKSolutionReader(const std::string& filename) : filename_(filename) {}

  CartesianDomainDefinition getCartesianDomain();

  void setSolution(const IndexArray& interior_indices, SolutionState& state);

 private:
  const std::string filename_;

  int dimension_;
};

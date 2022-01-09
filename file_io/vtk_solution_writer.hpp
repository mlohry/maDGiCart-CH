#pragma once

#include <string>

class SolutionState;
class Discretization2DCart;
class Discretization3DCart;
class SpatialDiscretization;


void write_solution_to_vtk(const SolutionState& state, const SpatialDiscretization& geom, std::string filename_withoutext);
void write_solution_to_vtk_2d(const SolutionState& state, const Discretization2DCart& geom, std::string filename_withoutext);
void write_solution_to_vtk_3d(const SolutionState& state, const Discretization3DCart& geom, std::string filename_withoutext);

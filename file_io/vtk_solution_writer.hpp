#pragma once

#include <string>

class SolutionState;
class Discretization2DCart;


void write_solution_to_vtk(const SolutionState& state, const Discretization2DCart& geom, std::string filename_withoutext);

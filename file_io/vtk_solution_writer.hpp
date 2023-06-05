#pragma once

#include <string>
#include <vector>


class SolutionState;
class Discretization2DCart;
class Discretization3DCart;
class SpatialDiscretization;
class TimeIntegrator;


void write_solution_to_vtk(
    const std::vector<std::pair<std::string, double>>& iteration_status,
    const SolutionState&                               state,
    const SpatialDiscretization&                       geom,
    std::string                                        filename_withoutext);
void write_solution_to_vtk_2d(
    const std::vector<std::pair<std::string, double>>& iteration_status,
    const SolutionState&                               state,
    const Discretization2DCart&                        geom,
    std::string                                        filename_withoutext);
void write_solution_to_vtk_3d(
    const std::vector<std::pair<std::string, double>>& iteration_status,
    const SolutionState&                               state,
    const Discretization3DCart&                        geom,
    std::string                                        filename_withoutext);

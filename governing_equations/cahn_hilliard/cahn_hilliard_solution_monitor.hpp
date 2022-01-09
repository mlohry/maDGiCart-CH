#pragma once

#include <string>
#include <utility>
#include <vector>

class TimeIntegrableRHS;
class SolutionState;
class Discretization2DCart;
class Discretization3DCart;

std::vector<std::pair<std::string, double>> cahn_hilliard_solution_monitor(
    const TimeIntegrableRHS& rhs,
    const Discretization2DCart& geom,
    const SolutionState&     state);

std::vector<std::pair<std::string, double>> cahn_hilliard_solution_monitor(
    const TimeIntegrableRHS& rhs,
    const Discretization3DCart& geom,
    const SolutionState&     state);

#pragma once

#include <string>
#include <utility>
#include <vector>

class TimeIntegrableRHS;
class SolutionState;
class Discretization2DCart;

std::vector<std::pair<std::string, double>> cahn_hilliard_solution_monitor(
    const TimeIntegrableRHS& rhs,
    const Discretization2DCart& geom,
    const SolutionState&     state);

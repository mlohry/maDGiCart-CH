#pragma once

#include <string>
#include <utility>
#include <vector>

class TimeIntegrableRHS;
class SolutionState;

std::vector<std::pair<std::string, double>> residual_calculator(
    const TimeIntegrableRHS& rhs,
    const SolutionState&     residuals);

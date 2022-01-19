#pragma once

#include <string>
#include <utility>
#include <vector>

class SolutionState;


class SolutionMonitor {
 public:
  virtual std::vector<std::pair<std::string, double>> eval(const SolutionState& state, double time, int iter) = 0;
};

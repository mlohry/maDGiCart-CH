#include "residual_calculator.hpp"

#include "governing_equations/time_integrable_rhs.hpp"

namespace {
double abs_res = 0;
double max_res = 0;
double rel_res = 0;
}

std::vector<std::pair<std::string, double>>
residual_calculator(const TimeIntegrableRHS& rhs, const SolutionState& residuals)
{
  std::vector<std::pair<std::string, double>> name_value_pairs;
  for (const auto& eqn_name : rhs.equationNames()) {
    name_value_pairs.emplace_back(std::pair<std::string, double>{"res_" + eqn_name, 1});
  }

  auto idx = read_access(rhs.interiorIndices());

  for (int ieqn = 0; ieqn < rhs.nEquations(); ++ieqn) {

    ReduceSumReal squared_norm(0);
    auto array = read_access(residuals.getVec(ieqn));

    maDGForAll(ii, 0, idx.size(), {
      squared_norm += pow(array[idx[ii]], 2.0);
    });
    name_value_pairs[ieqn].second = sqrt(squared_norm.get() / idx.size());
  }

  abs_res = name_value_pairs[0].second;
  max_res = std::max(max_res, abs_res);
  rel_res = abs_res / max_res;
  name_value_pairs.emplace_back(std::pair<std::string, double>{"rel_res", rel_res});


  return name_value_pairs;
}

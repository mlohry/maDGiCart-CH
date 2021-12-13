#include "residual_calculator.hpp"

#include "governing_equations/time_integrable_rhs.hpp"


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
    auto array = read_access(residuals.getVec(ieqn).asArray());

    maDGForAll(ii, 0, idx.size(), {
      squared_norm += pow(array[idx[ii]], 2.0);
    });
    name_value_pairs[ieqn].second = sqrt(squared_norm.get() / idx.size());
  }

  return name_value_pairs;
}

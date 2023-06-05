#include "constant_initial_conditions.hpp"
#include "data_structures/solution_state.hpp"
#include "time_integrable_rhs.hpp"


void
ConstantInitialConditions::set(const TimeIntegrableRHS& rhs, SolutionState& state) const
{
  auto         idx = read_access(rhs.interiorIndices());
  auto         c   = write_access(state.getVec(0));
  const double val = value_;

  maDGForAll(i, 0, idx.size(), { c[idx[i]] = val; });
}

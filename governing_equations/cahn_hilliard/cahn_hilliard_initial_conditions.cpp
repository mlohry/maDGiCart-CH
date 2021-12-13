#include "cahn_hilliard_initial_conditions.hpp"
#include "cahn_hilliard_parameters.hpp"
#include "governing_equations/time_integrable_rhs.hpp"
#include "data_structures/solution_state.hpp"

#include <random>


CahnHilliardInitialConditions::CahnHilliardInitialConditions(CahnHilliardParameters& params)
    : m_(params.m()), min_(params.initialMin()), max_(params.initialMax())
{
}


void
CahnHilliardInitialConditions::set(const TimeIntegrableRHS& rhs, SolutionState& state) const
{
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  //  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::mt19937                            gen(1);  // Standard mersenne_twister_engine seeded with 1
  std::uniform_real_distribution<real_wp> dist(min_, max_);


  auto idx = read_access_host(rhs.interiorIndices());
  auto c   = write_access_host(state.getVec(0).asArray());

  /**
   * Would prefer to use a maDGForAllHost here but std::uniform_real_distribution lambda copy
   * has to be mutable.
   */
  for (int i = 0; i < idx.size(); ++i) {
    c[idx[i]] = dist(gen);
  }
}
#pragma once

#include "data_structures/managed_array_2d.hpp"
#include "data_structures/solution_state.hpp"


class CahnHilliardState : public SolutionState {
 public:
  using BaseType = ManagedArray2DNonOwning<real_wp>;

  CahnHilliardState(const ManagedArrayOwner& owner, const std::string& state_name, int ni, int nj, int nhalo)
      : SolutionState(owner, state_name, 1, (ni + 2 * nhalo) * (nj + 2 * nhalo)),
        ni_(ni),
        nj_(nj),
        state_array_(const_cast<ManagedArray<real_wp>&>(getVec(0)), ni_, nj_, nhalo)
  {
  }

  std::vector<std::string> equationNames() const override { return {"c"}; }

  const BaseType& c() const
  {
    return state_array_;
  }

  BaseType& c()
  {
    return state_array_;
  }

 private:
  const int ni_;
  const int nj_;

  BaseType state_array_;
};

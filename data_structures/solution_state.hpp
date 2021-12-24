#pragma once

#include "data_structures/array_set.hpp"
#include "managed_array.hpp"


template<typename T>
class SolutionStateBase : public T {
 public:
  SolutionStateBase(
      const ManagedArrayOwner&  owner,
      const std::string&        state_name,
      int_t                     nvecs,
      int_t                     size)
      : T(owner, state_name, nvecs, size)
  {
  }

  virtual std::vector<std::string> equationNames() const = 0;

};

class SolutionState : public SolutionStateBase<ArraySet<ScalarArray>> {
 public:
  using SolutionStateBase<ArraySet<ScalarArray>>::SolutionStateBase;
};

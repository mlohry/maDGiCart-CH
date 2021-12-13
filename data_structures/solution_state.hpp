#pragma once


#include "managed_array.hpp"
#include "data_structures/mesh_datatypes.hpp"


template<typename T>
class SolutionStateBase : public T {
 public:
  SolutionStateBase(
      const ManagedArrayOwner&  owner,
      const std::string&        state_name,
      int_t                     nvecs,
      const StridedMatrixShape& shape)
      : T(owner, state_name, nvecs, shape)
  {
  }

  virtual std::vector<std::string> equationNames() const = 0;

};

class SolutionState : public SolutionStateBase<NodalScalarSet> {
 public:
  using SolutionStateBase<NodalScalarSet>::SolutionStateBase;
};

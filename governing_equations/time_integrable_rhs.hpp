#pragma once

#include <limits>
#include <string>
#include <vector>

#include "utils/observer.hpp"
#include "utils/noncopyable.hpp"

#include "typedefs.hpp"
#include "data_structures/solution_state.hpp"


class TimeIntegrableRHS : public Observable, private NonCopyable
{
 public:
  virtual ~TimeIntegrableRHS() = default;

  virtual void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) = 0;


  virtual std::vector<std::string> equationNames() const {
    std::vector<std::string> names;
    for (int_t i = 0; i < nEquations(); ++i )
      names.push_back("u"+std::to_string(i));
    return names;
  }

  virtual int_t nEquations() const = 0;

  virtual int_t dofsPerEquation() const = 0;

  virtual std::unique_ptr<SolutionState> createSolutionState() const = 0;

  virtual const IndexArray& interiorIndices() const = 0;
};

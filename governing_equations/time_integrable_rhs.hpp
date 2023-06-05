#pragma once

#include <limits>
#include <string>
#include <vector>

#include "utils/noncopyable.hpp"
#include "utils/observer.hpp"

#include "data_structures/csr_matrix.hpp"
#include "data_structures/solution_state.hpp"
#include "typedefs.hpp"

class SpatialDiscretization;

class TimeIntegrableRHS : public Observable, private NonCopyable {
 public:
  virtual ~TimeIntegrableRHS() = default;

  virtual void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) = 0;


  virtual std::vector<std::string> equationNames() const
  {
    std::vector<std::string> names;
    for (int_t i = 0; i < nEquations(); ++i)
      names.push_back("u" + std::to_string(i));
    return names;
  }

  virtual int_t nEquations() const = 0;

  virtual int_t dofsPerEquation() const = 0;

  virtual std::unique_ptr<SolutionState> createSolutionState() const = 0;

  virtual std::unique_ptr<CSRMatrix> createSparseMatrix() const
  {
    Logger::get().FatalMessage("TimeIntegrableRHS::createSparseMatrix not implemented");
    return nullptr;
  }

  virtual void evalJacobian(const SolutionState& flovars, double time, CSRMatrix& J)
  {
    Logger::get().FatalMessage("TimeIntegrableRHS::createSparseMatrix not implemented");
  }

  virtual std::vector<int> nNonZerosPerRow() const
  {
    Logger::get().FatalMessage("TimeIntegrableRHS::createSparseMatrix not implemented");
    return std::vector<int>{};
  }

  virtual const IndexArray& interiorIndices() const = 0;

  virtual std::unique_ptr<TimeIntegrableRHS> clone(SpatialDiscretization&) const
  {
    Logger::get().FatalMessage("TimeIntegrableRHS::clone not implemented");
    return nullptr;
  }

  virtual const SpatialDiscretization* hasGeometry() const
  {
    Logger::get().FatalMessage("TimeIntegrableRHS::hasGeometry not implemented");
    return nullptr;
  }

  virtual double cflScale() const {
    Logger::get().FatalMessage("TimeIntegrableRHS::cflScale not implemented");
    return 0;
  }
};

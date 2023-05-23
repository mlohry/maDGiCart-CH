#pragma once

#include "data_structures/managed_array_owner.hpp"
#include "logger/logger.hpp"


enum class BCType { Periodic, Neumann, Dirichlet };

struct CartesianDomainDefinition {
  int    nx;
  int    ny;
  int    nz;
  double xbeg;
  double xend;
  double ybeg;
  double yend;
  double zbeg;
  double zend;
  BCType xbc;
  BCType ybc;
  BCType zbc;
  int    nhalo;
};

class SolutionState;

class SpatialDiscretization : public ManagedArrayOwner {
 public:
  SpatialDiscretization(std::string classname, const CartesianDomainDefinition& domain)
      : ManagedArrayOwner(classname), domain_(domain)
  {
  }
  virtual ~SpatialDiscretization() = default;

  const CartesianDomainDefinition& domain() const { return domain_; }

  virtual std::unique_ptr<SpatialDiscretization> createCoarsenedDiscretization() const
  {
    Logger::get().FatalMessage("createCoarsenedDiscretization not implemented.");
    return nullptr;
  }

  virtual void interpolateFineToCoarse(
      const SolutionState&         fine_state,
      const SpatialDiscretization& coarse_geom,
      SolutionState&               coarse_state) const
  {
    Logger::get().FatalMessage("interpolateFineToCoarse not implemented.");
  }

  virtual void interpolateCoarseToFine(
      const SolutionState&         coarse_state,
      const SpatialDiscretization& fine_geom,
      SolutionState&               fine_state) const
  {
    Logger::get().FatalMessage("interpolateCoarseToFine not implemented.");
  }

 private:
  const CartesianDomainDefinition domain_;
};

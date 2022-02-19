#pragma once

#include "data_structures/managed_array_owner.hpp"


enum class BCType { Periodic, Neumann };

struct CartesianDomainDefinition {
  int    nx;
  int    ny;
  int    nz;
  double xbeg;
  double xend;
  double ybeg;
  double zbeg;
  BCType xbc;
  BCType ybc;
  BCType zbc;
  int nhalo;
};

class SpatialDiscretization : public ManagedArrayOwner {
 public:
  SpatialDiscretization(std::string classname, const CartesianDomainDefinition& domain)
      : ManagedArrayOwner(classname), domain_(domain)
  {
  }
  virtual ~SpatialDiscretization() = default;

 protected:
  const CartesianDomainDefinition& domain() const { return domain_; }

 private:
  const CartesianDomainDefinition domain_;
};

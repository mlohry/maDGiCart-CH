#pragma once

#include "data_structures/managed_array_owner.hpp"


struct CartesianDomainDefinition {
  int nx;
  int ny;
  int nz;
  double xbeg;
  double xend;
  double ybeg;
  double zbeg;
};



class SpatialDiscretization : public ManagedArrayOwner {
 public:
  SpatialDiscretization(std::string classname)
      : ManagedArrayOwner(classname)
  {}
  virtual ~SpatialDiscretization() = default;
};

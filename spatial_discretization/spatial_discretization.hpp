#pragma once

#include "data_structures/managed_array_owner.hpp"


class SpatialDiscretization : public ManagedArrayOwner {
 public:
  SpatialDiscretization(std::string classname)
      : ManagedArrayOwner(classname)
  {}
  virtual ~SpatialDiscretization() = default;
};

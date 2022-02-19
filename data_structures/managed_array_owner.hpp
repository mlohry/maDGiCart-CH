#pragma once

#include "utils/noncopyable.hpp"
#include "utils/uniquely_named_object.hpp"


class ManagedArrayOwner : private NonCopyable {
 public:
  ManagedArrayOwner() {}

  explicit ManagedArrayOwner(const std::string& name) : name_(name) {}

  virtual ~ManagedArrayOwner() = default;

  const std::string& objectName() const { return name_; }

 private:
  const std::string name_ = "UnnamedManagedArrayOwner";
};

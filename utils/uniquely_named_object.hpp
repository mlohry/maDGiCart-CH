#pragma once

#include <string>

#include <boost/core/demangle.hpp>


template <class T>
class UniquelyNamedObject
{
 public:
  UniquelyNamedObject()
      : name_(boost::core::demangle(typeid(T).name()) + "_" + std::to_string(nextID_++))
  {
  }

  const std::string& objectName() const { return name_; }

 private:
  std::string name_;
  static int        nextID_;
};

template <class T>
int UniquelyNamedObject<T>::nextID_ = 0;

#pragma once

#include <boost/core/demangle.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <typeinfo>


/**
  * A singleton class to handle registration of factories (or any function).
  * Using Scott Meyers' singleton pattern but with templated type.
  *
  * Bases classes must define a CtorType for constructor signature.
  */
template <class ClassType>
class FactoryRegistry
{
 public:
  using Factory = typename ClassType::CtorType;

  // Return the (singleton) Registry object
  static FactoryRegistry<ClassType>& get()
  {
    static FactoryRegistry<ClassType> instance;
    return instance;
  }

  // Declare a class factory by name
  bool add(const std::string &name, Factory factory)
  {
    registry.insert(std::make_pair(name, factory));
    return true;
  }

  // return the named factory function
  Factory lookup(const std::string& name) const
  {
    auto i = registry.find(name);
    if (i == registry.end()) {
      throw std::runtime_error("Unable to find " + name);
    }
    return i->second;
  }

  friend std::ostream& operator<<(std::ostream& os, const FactoryRegistry<ClassType>& r)
  {
    os << "Registry contents for type " << boost::core::demangle(typeid(ClassType).name()) << std::endl;
    for (const auto& i : r.registry){
      os << i.first << std::endl;
    }
    return os;
  }

 private:
  FactoryRegistry() {}
  ~FactoryRegistry() {}

  std::map<std::string, Factory> registry;
};

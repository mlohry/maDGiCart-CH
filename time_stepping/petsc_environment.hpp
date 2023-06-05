#pragma once

#include <string>

class PetscEnvironment
{
 public:
  PetscEnvironment(const std::string& petsc_cmd_line);
  ~PetscEnvironment();
};

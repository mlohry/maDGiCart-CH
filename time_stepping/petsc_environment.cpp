#include "petsc_environment.hpp"
#include "logger/logger.hpp"

#include <petscsys.h>
#include <iterator>
#include <string>
#include <vector>

namespace {

class PetscEnvironmentSingleton {
 public:
  static PetscEnvironmentSingleton& get(char** *argv, int *argc)
  {
    static PetscEnvironmentSingleton s(argv, argc);
    return s;
  }

  PetscEnvironmentSingleton(const PetscEnvironmentSingleton&) = delete;
  PetscEnvironmentSingleton& operator=(const PetscEnvironmentSingleton&) = delete;

 private:
  PetscEnvironmentSingleton(char*** argv, int *argc) {
    PetscInitialize(argc, argv, 0, 0);
  }
  ~PetscEnvironmentSingleton() {
    PetscFinalize();
  }
};

}


PetscEnvironment::PetscEnvironment(const std::string& petsc_cmd_line)
{
  std::istringstream iss(petsc_cmd_line);

  std::vector<std::string> petsc_args{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};

  petsc_args.insert(petsc_args.begin(), "maDG");  // executable name

  char** argv = new char*[petsc_args.size()];
  for (unsigned i = 0; i < petsc_args.size(); ++i) {
    argv[i] = const_cast<char*>(petsc_args[i].c_str());
  }
  int argc = petsc_args.size();


  std::string fullcmd;
  for (const auto& c : petsc_args) {
    fullcmd += c + " ";
  }

  // This also serves to ensure that Parallel has been initialized
  // so that petsc isn't initializing or destroying the MPI environment.
  Logger::get().InfoMessage("Initializing Petsc with command line options: " + fullcmd);
  PetscEnvironmentSingleton::get(&argv, &argc);
  Logger::get().InfoMessage("Petsc initialization complete.");

  delete[] argv;

  /// enables logging for PetscMemoryGetMaximumUsage
//  PetscMemorySetGetMaximumUsage();
}


PetscEnvironment::~PetscEnvironment()
{
  Logger::get().InfoMessage("Finalizing Petsc.");

}

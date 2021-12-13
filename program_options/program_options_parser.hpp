#pragma once

#include <boost/program_options.hpp>


class ProgramOptionsParser
{
 public:
  /*
   * Using a singleton here in combination with boost::program_options
   * creates a problem where re-parsed options throw errors. This is only
   * an issue in the automated tests where multiple solver instances get
   * initiated. Solution: don't use a singleton.
   */
  ProgramOptionsParser();
  ~ProgramOptionsParser() {}

  void parseInputOptions(const std::vector<std::string>& cmd_line);
  std::string getConfigurationFileTemplate() const;

 private:

  boost::program_options::options_description base{"Base options"};
  boost::program_options::options_description file{"File options"};
  boost::program_options::options_description physics{"Physics options"};
  boost::program_options::options_description discretization{"Spatial discretization options"};
  boost::program_options::options_description time_stepping{"Time stepping options"};
  boost::program_options::options_description visible;

  boost::program_options::variables_map vm;

  void printRegistryContents(std::ostream& out) const;
};

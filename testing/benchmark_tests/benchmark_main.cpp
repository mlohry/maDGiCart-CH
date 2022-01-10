#include "celero/Celero.h"
#include "logger/logger.hpp"

///
/// This is the main(int argc, char** argv) for the entire celero program.
/// You can write your own, or use this macro to insert the standard one into the project.
///

int
main(int argc, char** argv)
{
  Logger::get().setLogLevel(LogLevel::warning);
  celero::Run(argc, argv);
  return 0;
}

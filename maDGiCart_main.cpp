#include <string>
#include <vector>
#include "initialization/puppeteer.hpp"

int main(int argc, char* argv[])
{
  const std::vector<std::string> cmd_line(argv + 1, argv + argc);

  Puppeteer puppeteer(cmd_line);

  puppeteer.run();

  return 0;
}

#pragma once

#include <string>

class Machine {
 public:
 public:
  static Machine& get()
  {
    static Machine instance;
    return instance;
  }


  std::string getHostname() const;
  std::string getUsername() const;
  std::string getProcessorModel() const;
  std::string getDeviceModel() const;
  int         getDevicesPerNode() const;
  std::string getOpenMPInfo() const;

 private:
  Machine();
  ~Machine() = default;
};

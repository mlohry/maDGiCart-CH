#pragma once

#include <chrono>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>

#include "typedefs.hpp"


enum class LogLevel { trace, debug, info, warning, error, fatal };


class Logger {
 public:
  static Logger& get()
  {
    static Logger S;
    return S;
  }

  void disable();
  void enable();

  void               initLogFile(std::string file, LogLevel lvl = LogLevel::trace);
  const std::string& logFilename() const { return log_filename_; }
  void               closeLogFile();
  void               setLogLevel(LogLevel lvl);

  void updateLog();


  void setTimeMonitor(const std::vector<std::pair<std::string, double>>& name_value_pairs);
  void setResidualMonitor(const std::vector<std::pair<std::string, double>>& name_value_pairs);
  void setSolutionMonitor(const std::vector<std::pair<std::string, double>>& name_value_pairs);

  const std::vector<std::pair<std::string, double>> getTimeMonitor() const { return time_status_; }
  const std::vector<std::pair<std::string, double>> getResidualMonitor() const { return residual_status_; }
  const std::vector<std::pair<std::string, double>> getSolutionMonitor() const { return solution_status_; }

  void Message(const std::string& str, LogLevel lvl);
  void TraceMessage(std::string str);
  void DebugMessage(std::string str);
  void InfoMessage(std::string str);
  void WarningMessage(std::string str);
  void ErrorMessage(std::string str);
  void FatalMessage(std::string str);
  void WarningAssert(bool, const std::string&);
  void FatalAssert(bool, const std::string&);

  void logDeviceMemoryTransfer(const std::string& description);
  void closeDeviceMemoryTransferLog();


  class Timer {
   public:
    Timer(const std::string& msg, const LogLevel log_level);
    virtual ~Timer();
    real_t elapsed();

   protected:
    const std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    const std::string                                                 msg_;
    const LogLevel                                                    log_level_;
    bool                                                              elapsed_called_;
  };

  Timer timer(const std::string& msg, const LogLevel log_level = LogLevel::trace) { return Timer(msg, log_level); }


 private:
  Logger();
  ~Logger() = default;
  static std::string makeBanner();

  const std::string banner_;
  std::string       log_filename_;

  int_t n_update_calls_;

  bool initialized_;
  bool auto_flush_;
  bool log_level_overridden_ = false;

  std::unique_ptr<std::ofstream> memlog_;
  int_t                          n_device_to_host_copies_ = 0;

  std::vector<std::pair<std::string, double>> time_status_;
  std::vector<std::pair<std::string, double>> residual_status_;
  std::vector<std::pair<std::string, double>> solution_status_;

  boost::shared_ptr<boost::log::sinks::synchronous_sink<boost::log::sinks::text_file_backend>>    g_file_sink_;
  boost::shared_ptr<boost::log::sinks::synchronous_sink<boost::log::sinks::text_ostream_backend>> g_console_sink_;
};

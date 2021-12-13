#include "logger.hpp"
#include "git_info.hpp"

#include "parallel/machine.hpp"

#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/async_frontend.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/stacktrace.hpp>
#include <iomanip>
#include <iostream>
#include <tabulate/table.hpp>


namespace logging  = boost::log;
namespace src      = boost::log::sources;
namespace sinks    = boost::log::sinks;
namespace keywords = boost::log::keywords;


static const int colwidth = 13;

static tabulate::Table residual_header_table;


Logger::Logger()
    : banner_(makeBanner()),
      log_filename_("default"),
      n_update_calls_(0)
{
  if (!initialized_) {
    logging::add_common_attributes();  // enables timestamps and such

    // if the pointer is not null, a log file has already been initialized.
    // remove that sink pointer and initialize a new one..
    if (g_file_sink_) {
      g_file_sink_->flush();
      logging::core::get()->remove_sink(g_file_sink_);
      g_file_sink_.reset();
    }

    if (g_console_sink_) {
      g_console_sink_->flush();
      logging::core::get()->remove_sink(g_console_sink_);
      g_console_sink_.reset();
    }

    g_console_sink_ = logging::add_console_log(std::clog, keywords::format = "%Message%");

    initialized_ = true;
  }
}


void
Logger::setLogLevel(LogLevel lvl)
{
  switch (lvl) {
    case LogLevel::trace:
      logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::trace);
      break;
    case LogLevel::debug:
      logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::debug);
      break;
    case LogLevel::info:
      logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);
      break;
    case LogLevel::warning:
      logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::warning);
      break;
    case LogLevel::error:
      logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::error);
      break;
    case LogLevel::fatal:
      logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::fatal);
      break;
  }

  log_level_overridden_ = true;
}

void
Logger::disable()
{
  logging::core::get()->set_logging_enabled(false);
}

void
Logger::enable()
{
  logging::core::get()->set_logging_enabled(true);
}

void
Logger::initLogFile(std::string file, LogLevel lvl)
{
  log_filename_ = file;

  // if the pointer is not null, a log file has already been initialized.
  // remove that sink pointer and initialize a new one..
  if (g_file_sink_) {
    g_file_sink_->flush();
    logging::core::get()->remove_sink(g_file_sink_);
    g_file_sink_.reset();
  }

  g_file_sink_ = logging::add_file_log(
      keywords::target        = "",  // Forces logger to not overwrite a log file, increment the file number instead.
      keywords::file_name     = log_filename_ + "_%N.log",
      keywords::rotation_size = 1024 * 1024 * 500,  // 500MB log file max
      keywords::auto_flush    = false,
      keywords::format        = "%Message% #[%TimeStamp%]");


  if (!log_level_overridden_) {
    setLogLevel(lvl);
  }
  BOOST_LOG_TRIVIAL(info) << banner_;
}

void
Logger::updateLog()
{

  if (n_update_calls_ == 0) {

    std::vector<std::string> logger_output_strings;

    for (const auto& timepair : time_status_){
      logger_output_strings.push_back(timepair.first);
    }

    for (const auto& residpair : residual_status_){
      logger_output_strings.push_back(residpair.first);
    }

    for (const auto& solpair : solution_status_){
      logger_output_strings.push_back(solpair.first);
    }


    BOOST_LOG_TRIVIAL(warning) << "#\n# Logger outputs: ";

    std::string outputs;
    for (std::size_t i = 0; i < logger_output_strings.size(); ++i) {
      outputs += logger_output_strings[i] + " ";
    }
    BOOST_LOG_TRIVIAL(warning) << outputs;

    residual_header_table = tabulate::Table();
    residual_header_table.add_row({logger_output_strings.begin(), logger_output_strings.end()});

    residual_header_table.format().border_top(" ").border_right(" ").border_bottom(" ").border_left(" ");

    for (auto& cell : residual_header_table[0]) {
      cell.format()
          .width(colwidth - 1)
          .font_align(tabulate::FontAlign::right)
          .border_left(" ")
          .corner_top_left(" ")
          .corner_bottom_left(" ")
          .corner_top_right(" ")
          .corner_bottom_right(" ");
    }

    residual_header_table[0][0].format().border_left("#").corner_top_left("#").corner_bottom_left("#").width(
        colwidth - 1);
  }

  if (n_update_calls_ == 0 || n_update_calls_ % 50 == 0) {
    BOOST_LOG_TRIVIAL(warning) << residual_header_table;
  }

  std::ostringstream ssres;
  ssres << " ";
  for (const auto& timepair : time_status_){
    ssres << std::scientific << std::setprecision(5) << std::setw(colwidth) << timepair.second;
  }

  for (const auto& residpair : residual_status_){
    ssres << std::scientific << std::setprecision(5) << std::setw(colwidth) << residpair.second;
  }

  for (const auto& solpair : solution_status_){
    ssres << std::scientific << std::setprecision(5) << std::setw(colwidth) << solpair.second;
  }

  BOOST_LOG_TRIVIAL(warning) << ssres.str();
  if (g_file_sink_ && auto_flush_) {
    g_file_sink_->flush();
  }
  n_update_calls_++;
}


void
Logger::FatalAssert(bool condition, const std::string& str)
{
  if (!condition) {
    FatalMessage(str);
  }
}

void
Logger::WarningAssert(bool condition, const std::string& str)
{
  if (!condition) {
    WarningMessage(str);
  }
}

void
Logger::Message(const std::string& str, LogLevel lvl)
{
  switch (lvl) {
    case LogLevel::trace:
      TraceMessage(str);
      break;
    case LogLevel::debug:
      DebugMessage(str);
      break;
    case LogLevel::info:
      InfoMessage(str);
      break;
    case LogLevel::warning:
      WarningMessage(str);
      break;
    case LogLevel::error:
      ErrorMessage(str);
      break;
    case LogLevel::fatal:
      FatalMessage(str);
      break;
  }
}

void
Logger::FatalMessage(std::string str)
{
  // writing to standard output here to catch logged messages from non-master processes.
  // other logger messages are usually called by all processes.
  std::cout << "# [FATAL] " << str << std::endl;

  str = "**************************************\n\n" + str + "\n\n**************************************";
  boost::replace_all(str, "\n", "\n# [FATAL] ");
  BOOST_LOG_TRIVIAL(fatal) << "# [FATAL] " << str;
  if (g_file_sink_) {
    g_file_sink_->flush();
  }
  abort();
}

void
Logger::ErrorMessage(std::string str)
{
  boost::replace_all(str, "\n", "\n# [Trace] ");  // replaces newlines to newlines with comment character
  BOOST_LOG_TRIVIAL(error) << "# [Trace] " << str;
}

void
Logger::TraceMessage(std::string str)
{
  boost::replace_all(str, "\n", "\n# [Trace] ");  // replaces newlines to newlines with comment character
  BOOST_LOG_TRIVIAL(trace) << "# [Trace] " << str;
}

void
Logger::DebugMessage(std::string str)
{
  boost::replace_all(str, "\n", "\n# [Debug] ");  // replaces newlines to newlines with comment character
  BOOST_LOG_TRIVIAL(debug) << "# [Debug] " << str;
}

void
Logger::WarningMessage(std::string str)
{
  boost::replace_all(str, "\n", "\n# [Warning] ");  // replaces newlines to newlines with comment character
  BOOST_LOG_TRIVIAL(warning) << "# [Warning] " << str;
}


void
Logger::InfoMessage(std::string str)
{
  boost::replace_all(str, "\n", "\n# [Info] ");  // replaces newlines to newlines with comment character
  BOOST_LOG_TRIVIAL(info) << "# [Info] " << str;
}



void
Logger::closeLogFile()
{
  if (g_file_sink_) {
    g_file_sink_->flush();
    logging::core::get()->remove_sink(g_file_sink_);
    g_file_sink_.reset();
  }

  if (g_console_sink_) {
    g_console_sink_->flush();
    logging::core::get()->remove_sink(g_console_sink_);
    g_console_sink_.reset();
  }

  initialized_         = false;
  n_update_calls_      = 0;
}


void
Logger::setTimeMonitor(const std::vector<std::pair<std::string, double>>& name_value_pairs)
{
  time_status_ = name_value_pairs;
}

void
Logger::setSolutionMonitor(const std::vector<std::pair<std::string, double>>& name_value_pairs)
{
  solution_status_ = name_value_pairs;
}

void
Logger::setResidualMonitor(const std::vector<std::pair<std::string, double>>& name_value_pairs)
{
  residual_status_ = name_value_pairs;
}


Logger::Timer::Timer(const std::string& msg, const LogLevel log_level)
    : start_(std::chrono::high_resolution_clock::now()), msg_(msg), log_level_(log_level), elapsed_called_(false)
{
  Logger::get().Message("Timing: " + msg_ + "...", log_level_);
}

Logger::Timer::~Timer()
{
  if (!elapsed_called_) {
    elapsed();
  }
}

real_t
Logger::Timer::elapsed()
{
  std::chrono::duration<real_t> dt   = std::chrono::high_resolution_clock::now() - start_;
  real_t                        time = (real_t)dt.count();
  Logger::get().Message("Timing: " + msg_ + " finished in " + std::to_string(time) + " seconds.", log_level_);
  elapsed_called_ = true;
  return time;
}


void
Logger::logDeviceMemoryTransfer(const std::string& description)
{
  n_device_to_host_copies_++;
  if (!memlog_) {
    memlog_ = std::make_unique<std::ofstream>(logFilename() + "_memory.log");
  }
  *memlog_ << description << std::endl << boost::stacktrace::stacktrace() << std::endl;
}


void
Logger::closeDeviceMemoryTransferLog()
{
  if (memlog_) {
    *memlog_ << "Device-to-host copies: " + std::to_string(n_device_to_host_copies_) << std::endl;
    memlog_->close();
    memlog_.reset(nullptr);
  }

  InfoMessage("Device-to-host copies: " + std::to_string(n_device_to_host_copies_));
}


std::string
Logger::makeBanner()
{
  // clang-format off
  const std::string banner =
      R"(#                                           )""\n"
      R"(#                       ______ _____        )""\n"
      R"(#                       |  _  \  __ \       )""\n"
      R"(#        _ __ ___   __ _| | | | |  \/       )""\n"
      R"(#       | '_ ` _ \ / _` | | | | | __        )""\n"
      R"(#       | | | | | | (_| | |/ /| |_\ \       )""\n"
      R"(#       |_| |_| |_|\__,_|___/  \____/       )""\n"
      R"(#                                           )""\n";

  boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();

  std::stringstream ss;
  ss << banner;
  ss << R"(#     Hostname           )" << Machine::get().getHostname() << "\n";
  ss << R"(#     Username           )" << Machine::get().getUsername() << "\n";
  ss << R"(#     Work directory     )" << boost::filesystem::current_path().string() << "\n";
  ss << R"(#     Processor model    )" << Machine::get().getProcessorModel() << "\n";
  ss << R"(#     GPU model          )" << Machine::get().getDeviceModel() << "\n";
  ss << R"(#     GPUs per node      )" << Machine::get().getDevicesPerNode() << "\n";
  ss << R"(#     OpenMP max threads )" << Machine::get().getOpenMPInfo() << "\n";
  ss << R"(#     System time        )" << timeLocal << "\n";
  ss << R"(#     git branch         )" << GIT_BRANCH << "\n";
  ss << R"(#     git commit hash    )" << GIT_COMMIT_HASH << "\n";
  ss << R"(#     C++ compiler       )" << CXX_COMPILER_STRING << "\n";
  ss << R"(#                        )" << "\n";
  ss << R"(#)";

  // clang-format on

  return ss.str();
}

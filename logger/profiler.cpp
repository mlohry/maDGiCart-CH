#include "profiler.hpp"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <tabulate/table.hpp>


#include "logger.hpp"


std::string Profiler::finalize()
{
  tabulate::Table table;
  table.add_row({"Function", "# calls", "Time (s)", "Percentage"});

  struct Report {
    std::string name;
    long long   ncalls;
    double      time;
    double      percentage;
  };

  std::vector<Report> reports;
  double              max_time = 0;
  for (const auto& profile : function_profiles_) {
    max_time = std::max(max_time, profile.second.total_time_);
  }

  for (const auto& profile : function_profiles_) {
    reports.push_back(
        {profile.first,
         profile.second.ncalls_,
         profile.second.total_time_,
         profile.second.total_time_ / max_time * 100.0});
  }

  std::sort(reports.begin(), reports.end(), [](const Report& a, const Report& b) {
    return a.percentage > b.percentage;
  });


  auto fmt_float = [](real_t r) {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(3) << r;
    return ss.str();
  };

  auto fmt_pct = [](real_t r) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << r << "%";
    return ss.str();
  };

  for (const auto& report : reports) {
    table.add_row(
        {report.name,
         std::to_string(report.ncalls),
         fmt_float(report.time),
         fmt_pct(report.percentage)});
  }

  table.column(0).format().width(70);
  table.column(3).format().font_align(tabulate::FontAlign::right);

  std::stringstream ss;
  ss << table << std::endl;
  return ss.str();
}


void Profiler::tic(const std::string& function_name)
{
  function_profiles_[function_name].ncalls_++;
}

void Profiler::toc(const std::string& function_name, double elapsed)
{
  function_profiles_[function_name].total_time_ += elapsed;
}


Profiler::FunctionTimer::FunctionTimer(Profiler& profiler, const std::string& function_name)
    : profiler_(profiler),
      start_(std::chrono::high_resolution_clock::now()),
      function_name_(function_name)
{
  profiler_.tic(function_name_);
}

Profiler::FunctionTimer::~FunctionTimer()
{
  std::chrono::duration<double> dt      = std::chrono::high_resolution_clock::now() - start_;
  const double                  elapsed = (double)dt.count();
  profiler_.toc(function_name_, elapsed);
}

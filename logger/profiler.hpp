#pragma once


#include <boost/current_function.hpp>
#include <chrono>
#include <string>
#include <unordered_map>
#include "typedefs.hpp"


#define profile() \
  auto function_timer = Profiler::FunctionTimer(Profiler::get(), BOOST_CURRENT_FUNCTION);

#define profile_tic(UNIQUE_NAME) \
  auto function_timer_##UNIQUE_NAME = \
 std::make_unique<Profiler::FunctionTimer>(Profiler::get(), #UNIQUE_NAME);

#define profile_toc(UNIQUE_NAME) \
  function_timer_##UNIQUE_NAME.reset(nullptr);

class Profiler
{
 public:
  static Profiler& get()
  {
    static Profiler S;
    return S;
  }


  class FunctionTimer
  {
   public:
    FunctionTimer(Profiler& profiler, const std::string& function_name);
    ~FunctionTimer();

   private:
    Profiler&                                                         profiler_;
    const std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    const std::string                                                 function_name_;
  };

  std::string finalize();

 private:
  Profiler() = default;
  ~Profiler() = default;


  void tic(const std::string& function_name);
  void toc(const std::string& function_name, double elapsed);


  struct FunctionProfile {
    long long ncalls_     = 0;
    double    total_time_ = 0;
  };

  std::unordered_map<std::string, FunctionProfile> function_profiles_;
};

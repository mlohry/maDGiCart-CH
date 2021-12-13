#pragma once

#include <functional>
#include <map>
#include <vector>

#include "logger/profiler.hpp"

/// event classifications mapping to the observer pattern
enum class Event {
  TimeStart,
  TimeComplete,
  TimeStepStart,
  TimeStepComplete,
  InnerIterationStart,
  InnerIterationComplete,
  NonlinearIterationStart,
  NonlinearIterationComplete,
  NonlinearSolveCompletion,
  LinearStepCompletion,
  LinearSolveCompletion,
  RHSEvaluation
};


/**
 * Mix-in class to make a class "observable" and
 * allow the registration of observers.
 *
 * Modified from
 * https://juanchopanzacpp.wordpress.com/2013/02/24/simple-observer-pattern-implementation-c11/
 */
class Observable
{
 public:
  virtual ~Observable() {}

  template <typename Observer>
  void registerObserver(const Event& event, Observer&& observer)
  {
    observers_[event].push_back(std::forward<Observer>(observer));
  }

  template <typename Observer>
  void registerObserverFront(const Event& event, Observer&& observer)
  {
    observers_[event].insert(observers_[event].begin(), std::forward<Observer>(observer));
  }

  /// notify specific observers corresponding to the event
  void notifyObservers(const Event& event) const
  {
    profile();
    /*
     * Note: for (const auto& obs : observers_.at(event))
     * will throw an exception if that event doesn't exist in the map.
     */
    auto it = observers_.find(event);
    if (it != observers_.end())
    {
      for (const auto& obs : it->second)
      {
        obs();
      }
    }
  }

  void copyObservers(const Observable& source)
  {
    for (const auto& kv : source.observers_) {
      for (const auto& f : kv.second) {
        this->observers_[kv.first].push_back(f);
      }
    }
  }

 private:
  std::map<Event, std::vector<std::function<void()>>> observers_;
};

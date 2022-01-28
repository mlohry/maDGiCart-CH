#pragma once

#include <functional>
#include <map>
#include <vector>

#include "logger/profiler.hpp"

/// event classifications mapping to the observer pattern
enum class Event { TimeStepComplete, SolutionUpdate };


/**
 * Mix-in class to make a class "observable" and
 * allow the registration of observers.
 *
 * Modified from
 * https://juanchopanzacpp.wordpress.com/2013/02/24/simple-observer-pattern-implementation-c11/
 */
class Observable {
 public:
  Observable()
  {
    observer_update_freq_[Event::TimeStepComplete]     = 1;
    observer_update_freq_[Event::SolutionUpdate]       = 1;
    observer_n_notifications_[Event::TimeStepComplete] = 0;
    observer_n_notifications_[Event::SolutionUpdate]   = 0;
  }

  virtual ~Observable() {}

  template <typename Observer>
  void registerObserver(const Event& event, Observer&& observer)
  {
    observers_[event].push_back(std::forward<Observer>(observer));
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
    if (it != observers_.end()) {
      if (observer_n_notifications_.at(event) % observer_update_freq_.at(event) == 0) {
        for (const auto& obs : it->second) {
          obs();
        }
      }
      observer_n_notifications_.at(event)++;
    }
  }

  void setNotifyFrequency(const Event& event, int freq) { observer_update_freq_[event] = freq; }

 private:
  std::map<Event, std::vector<std::function<void()>>> observers_;
  std::map<Event, int>                                observer_update_freq_;
  mutable std::map<Event, int>                        observer_n_notifications_;
};

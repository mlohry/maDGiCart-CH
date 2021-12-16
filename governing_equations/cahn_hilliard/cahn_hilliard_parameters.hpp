#pragma once

#include "program_options/program_options.hpp"

class CahnHilliardParameters {
 public:
  double m() const { return m_; }
  double sigma() const { return sigma_; }
  double eps2() const { return eps2_; }
  double initialMin() const { return initial_value_min_; }
  double initialMax() const { return initial_value_max_; }

 private:
  const double m_                 = Options::get().ch_m();
  const double sigma_             = Options::get().ch_sigma();
  const double eps2_              = Options::get().ch_eps2();
  const double initial_value_min_ = -0.005;
  const double initial_value_max_ = 0.005;
};
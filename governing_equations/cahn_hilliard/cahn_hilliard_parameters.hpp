#pragma once

class CahnHilliardParameters {
 public:
  double m() const { return m_; }
  double sigma() const { return sigma_; }
  double eps2() const { return eps2_; }
  double initialMin() const { return initial_value_min_; }
  double initialMax() const { return initial_value_max_; }

 private:
  const double m_                 = 0;
  const double sigma_             = 291.4488109293305;
  const double eps2_              = 2.419753086419753e-05;
  const double initial_value_min_ = -0.005;
  const double initial_value_max_ = 0.005;
};
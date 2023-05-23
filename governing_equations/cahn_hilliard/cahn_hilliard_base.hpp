#pragma once

#include "cahn_hilliard_parameters.hpp"
#include "governing_equations/time_integrable_rhs.hpp"


class CahnHilliardBase : public TimeIntegrableRHS {
 public:
  CahnHilliardBase(const CahnHilliardParameters& params)
      : params_(params),
        m_(params.m()),
        sigma_(params.sigma()),
        eps2_(params.eps2())
  {
  }

  int_t nEquations() const override { return 1; }

  std::vector<std::string> equationNames() const override { return {{"c"}}; }

 protected:
  double m() const { return m_; }
  double sigma() const { return sigma_; }
  double eps2() const { return eps2_; }

  const CahnHilliardParameters& params() const { return params_; }

 private:
  const CahnHilliardParameters params_;
  const double m_;
  const double sigma_;
  const double eps2_;
};

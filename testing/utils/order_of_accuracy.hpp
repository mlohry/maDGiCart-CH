#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <iomanip>

#include "typedefs.hpp"


inline double errorRateByNdof(double ndof0, double ndof1, double e0, double e1)
{
  return -log10(e1 / e0) / log10(ndof1 / ndof0);
}


class OrderOfAccuracy
{
 public:
  using dpair = std::pair<double, double>;

  void addSolutionError(std::size_t order, double ndof, double error)
  {
    solutions_[order].push_back(std::make_pair(ndof, error));

    std::sort(
        solutions_[order].begin(), solutions_[order].end(), [](const dpair& p0, const dpair& p1) {
          return p0.first < p1.first;
        });
  }


  std::map<std::size_t, std::vector<double>> getConvergenceRate() const
  {
    std::map<std::size_t, std::vector<double>> rate;

    for (const auto& kv : solutions_) {
      rate[kv.first].push_back(0);
      for (std::size_t i = 0; i < kv.second.size() - 1; ++i) {
        rate[kv.first].push_back(
            errorRateByNdof(
                kv.second[i].first,
                kv.second[i + 1].first,
                kv.second[i].second,
                kv.second[i + 1].second));
      }
    }

    return rate;
  }

  const std::map<std::size_t, std::vector<dpair>>& get() const { return solutions_; }


  real_t getMinimumError() const
  {
    real_t err = std::numeric_limits<real_t>::max();
    for (const auto& kv : solutions_) {
      for (std::size_t i = 0; i < kv.second.size(); ++i) {
        err = std::min(err, kv.second[i].second);
      }
    }
    return err;
  }

  void writeToFiles(const std::string& prefix)
  {
    for (const auto& kv : solutions_) {
      const std::string filename = prefix + "_N" + std::to_string(kv.first) + ".dat";
      std::ofstream     fout;
      fout.open(filename);

      for (std::size_t i = 0; i < kv.second.size(); ++i) {
        fout << std::scientific << std::setprecision(6) << kv.second[i].first << " "
             << kv.second[i].second << "\n";
      }

      fout.close();
    }
  }

 private:
  std::map<std::size_t, std::vector<dpair>> solutions_;
  friend std::ostream& operator<<(std::ostream&, const OrderOfAccuracy&);
};


inline std::ostream& operator<<(std::ostream& os, const OrderOfAccuracy& accuracy)
{
  os << std::setw(11) << "P-Order" << std::setw(11) << "nDOFs" << std::setw(11) << "error"
     << std::setw(11) << "order\n";
  os << std::string(44, '-') << "\n";

  for (const auto& kv : accuracy.solutions_) {
    // std::sort(kv.second.begin(), kv.second.end(),
    //           [](const dpair& p0, const dpair& p1){return p0.first < p1.first;});

    for (std::size_t i = 0; i < kv.second.size(); ++i) {
      if (i == 0) {
        os << std::setw(11) << kv.first;
        os << std::scientific << std::setprecision(3) << std::setw(11) << kv.second[i].first;
        os << std::scientific << std::setprecision(3) << std::setw(11) << kv.second[i].second
           << "\n";
      } else {
        double err = errorRateByNdof(
            kv.second[i - 1].first,
            kv.second[i].first,
            kv.second[i - 1].second,
            kv.second[i].second);
        os << std::setw(11) << kv.first;
        os << std::scientific << std::setprecision(3) << std::setw(11) << kv.second[i].first;
        os << std::scientific << std::setprecision(3) << std::setw(11) << kv.second[i].second;
        os << std::defaultfloat << std::setw(11) << err << "\n";
      }
    }
    os << std::string(44, '-') << "\n";
  }

  return os;
}

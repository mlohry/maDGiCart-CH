#include "asymptotic_limit.hpp"

#include <lsqcpp.h>

#include <iostream>
#include <vector>




static std::deque<double> time_history_;
static std::deque<double> val_history_;
static bool first_iter_ = true;
static Eigen::VectorXd last_solution_;

namespace {
std::vector<double> getLastN(const std::deque<double>& deq, int N)
{
  std::vector<double> vec;
  const size_t dsize = deq.size();
  for (size_t i = dsize-1-N; i < dsize; ++i){
    vec.push_back(deq[i]);
  }
  return vec;
}
}


class AsymptoteEstimatorP2Q2 {
 public:
  AsymptoteEstimatorP2Q2() {}

  Eigen::VectorXd initialGuess() const {
//    if (first_iter_) {
      //Eigen::VectorXd guess = Eigen::VectorXd::Constant(5,0.0);
      //guess(0) = val_history_.back();
      //return guess;

      Eigen::MatrixXd V(6,6);
      const auto x = getLastN(time_history_, 5);
      const auto y = getLastN(val_history_, 5);
      Eigen::VectorXd yvec(5);

      for (int i = 0; i < 5; ++i){
        V(i, 0) = 1.0;
        V(i, 1) = x[i];
        V(i, 2) = pow(x[i], 2.0);
        V(i, 3) = pow(x[i], 3.0);
        V(i, 4) = pow(x[i], 4.0);
        V(i, 5) = pow(x[i], 5.0);
        yvec(i) = y[i];
      }

      const Eigen::VectorXd coefs = (V.transpose()*V).inverse() * yvec;
      std::cout << coefs << "\n";

      return coefs;



//      return Eigen::VectorXd::Ones(5);
//    }
//    return last_solution_;
  }

  double asymptote(const Eigen::VectorXd& xval) const {
    first_iter_ = false;
    last_solution_ = xval;
    return xval(2) / xval(4);
  }

  void setValues(const std::deque<double>& timevec, const std::deque<double>& valvec)
  {
    time_history_ = timevec;
    val_history_ = valvec;
  }

  struct ErrorCaclulator {
    void operator()(const Eigen::VectorXd& xval, Eigen::VectorXd& fval, Eigen::MatrixXd& Jac) const
    {
      const double a0 = xval(0);
      const double a1 = xval(1);
      const double a2 = xval(2);
      const double b1 = xval(3);
      const double b2 = xval(4);

      fval.resize(val_history_.size());

      Jac.resize(val_history_.size() , 5);

      for (int i = 0; i < fval.size(); ++i){
        // d/dA0
        auto ddA0 = [&](double x){
          return 1.0 / (1.0+b1*x+b2*x*x);
        };
        auto ddA1 = [&](double x){
          return x / (1.0+b1*x+b2*x*x);
        };
        auto ddA2 = [&](double x){
          return x*x / (1.0+b1*x+b2*x*x);
        };
        auto ddB1 = [&](double x){
          return -x*(a0+x*(a1+a2*x))/pow(1.0+b1*x+b2*x*x, 2.0);
        };
        auto ddB2 = [&](double x){
          return -x*x*(a0+x*(a1+a2*x))/pow(1.0+b1*x+b2*x*x, 2.0);
        };

        const double x = time_history_[i];
        Jac(i, 0) = ddA0(x);
        Jac(i, 1) = ddA1(x);
        Jac(i, 2) = ddA2(x);
        Jac(i, 3) = ddB1(x);
        Jac(i, 4) = ddB2(x);
      }


      for (lsq::Index i = 0; i < fval.size(); ++i) {
        const double x      = time_history_[i];
        const double fitval = (a0 + a1 * x + a2 * x * x) / (1.0 + b1 * x + b2 * x * x);
        fval(i)             = fitval - val_history_[i];
      }
    }
  };

 private:
};


double
asymptoticLimit(const std::deque<double>& time, const std::deque<double>& value)
{
  using RationalType = AsymptoteEstimatorP2Q2;
  RationalType estimator;
  estimator.setValues(time, value);

  lsq::LevenbergMarquardt<double, RationalType::ErrorCaclulator> optimizer;
//  lsq::GradientDescent<double, RationalType::ErrorCaclulator> optimizer;
//  lsq::GaussNewton<double, RationalType::ErrorCaclulator, lsq::ArmijoBacktracking<double>> optimizer;

  optimizer.setMaxIterations(100);
  optimizer.setVerbosity(0);
  optimizer.setMinError(1.e-10);
//  optimizer.setMaxIterationsLM(100);
//  optimizer.setLambdaIncrease(2);
//  optimizer.setLambdaDecrease(0.5);
//  optimizer.setVerbosity(2);
//  optimizer.setThreads(0);

  auto result = optimizer.minimize(estimator.initialGuess());

//  std::cout << "Done! Converged: " << (result.converged ? "true" : "false") << " Iterations: " << result.iterations
//            << std::endl;

  return estimator.asymptote(result.xval);
}

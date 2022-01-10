#include "celero/Celero.h"
#include "data_structures/exec_includes.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_2d_fd.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"


class CahnHilliard2DFixture : public celero::TestFixture {
 public:
  CahnHilliard2DFixture() {}


  virtual std::vector<celero::TestFixture::ExperimentValue> getExperimentValues() const override
  {
    std::vector<celero::TestFixture::ExperimentValue> problemSpace;

    problemSpace.push_back({int64_t(64)});
    problemSpace.push_back({int64_t(128)});
    problemSpace.push_back({int64_t(256)});

    return problemSpace;
  }


  /// Before each run
  virtual void setUp(const celero::TestFixture::ExperimentValue& experimentValue)
  {
    this->gridsize_ = experimentValue.Value;
    geom_           = std::make_unique<Discretization2DCart>(this->gridsize_, 2, 0, 1, 0);

    CahnHilliardParameters ch_params;
    ch_rhs_ = std::make_unique<CahnHilliard2DFD>(*geom_, ch_params);

    state_     = ch_rhs_->createSolutionState();
    dstate_dt_ = ch_rhs_->createSolutionState();

    auto s = write_access_host(state_->getVec(0));
    maDGForAllHost(i, 0, s.size(), { s[i] = 0; });
  }


  int                                   gridsize_;
  std::unique_ptr<Discretization2DCart> geom_;
  std::unique_ptr<CahnHilliard2DFD>     ch_rhs_;
  std::unique_ptr<SolutionState>        state_;
  std::unique_ptr<SolutionState>        dstate_dt_;
};

BASELINE_F(CahnHilliard2D, CHRHS2D, CahnHilliard2DFixture, 10, 10)
{  //
  ch_rhs_->evalRHSImpl(*state_, 0, *dstate_dt_);
}

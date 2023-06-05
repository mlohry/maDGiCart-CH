#include "celero/Celero.h"
#include "data_structures/exec_includes.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_3d_fd.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"


class CahnHilliard3DFixture : public celero::TestFixture {
 public:
  CahnHilliard3DFixture() {}


  virtual std::vector<celero::TestFixture::ExperimentValue> getExperimentValues() const override
  {
    std::vector<celero::TestFixture::ExperimentValue> problemSpace;

    problemSpace.push_back({int64_t(64)});
    problemSpace.push_back({int64_t(128)});
    problemSpace.push_back({int64_t(256)});

    return problemSpace;
  }


  /// Before each run
  void setUp(const celero::TestFixture::ExperimentValue& experimentValue) override
  {
    this->gridsize_ = experimentValue.Value;

    CartesianDomainDefinition domain;
    domain.nx = domain.ny = domain.nz = this->gridsize_;
    domain.xbc = domain.ybc = domain.zbc = BCType::Periodic;

    domain.xbeg  = 0;
    domain.xend  = 1;
    domain.ybeg  = 0;
    domain.zbeg  = 0;
    domain.nhalo = 2;

    geom_ = std::make_unique<Discretization3DCart>(domain);

    CahnHilliardParameters ch_params;
    ch_rhs_ = std::make_unique<CahnHilliard3DFD>(*geom_, ch_params);

    state_     = ch_rhs_->createSolutionState();
    dstate_dt_ = ch_rhs_->createSolutionState();

    auto s = write_access_host(state_->getVec(0));
    maDGForAllHost(i, 0, s.size(), { s[i] = 0; });
  }


  int                                   gridsize_;
  std::unique_ptr<Discretization3DCart> geom_;
  std::unique_ptr<CahnHilliard3DFD>     ch_rhs_;
  std::unique_ptr<SolutionState>        state_;
  std::unique_ptr<SolutionState>        dstate_dt_;
};

BASELINE_F(CahnHilliard3D, CHRHS3D, CahnHilliard3DFixture, 0, 5)
{  //
  ch_rhs_->evalRHSImpl(*state_, 0, *dstate_dt_);
}

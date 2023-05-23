#include "cahn_hilliard_3d_fd.hpp"

#include "logger/logger.hpp"
#include "program_options/program_options.hpp"


CahnHilliard3DFD::CahnHilliard3DFD(Discretization3DCart& geom, const CahnHilliardParameters& params)
    : CahnHilliardBase(params), geom_(geom), kernel_variant_(Options::get().kernel_variant())
{
  laplacian_rhs_term_ = geom_.createRealArray();
  laplacian_argument_ = geom_.createRealArray();
  biharmonic_term_    = geom_.createRealArray();
  linear_term_        = geom_.createRealArray();

  Logger::get().FatalAssert(kernel_variant_ >= 0 && kernel_variant_ <= 2, "kernel_variant must be between 0 and 2.");
}


void
CahnHilliard3DFD::evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs)
{
//  auto timer = Logger::get().timer("CahnHilliard3DFD::evalRHSImpl");
  switch (kernel_variant_) {
    case 0:
      evalRHSBasic(flovars, time, rhs);
      break;
    case 1:
      evalRHSFused(flovars, time, rhs);
      break;
    case 2:
      evalRHSFullyFused(flovars, time, rhs);
      break;
  }
}


void
CahnHilliard3DFD::evalRHSBasic(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const ScalarSolutionState3D& state = dynamic_cast<const ScalarSolutionState3D&>(flovars);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));
  geom_.applyNeumannBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));

  /**
   * Compute  laplacian(u*c^3 -b*c)
   */
  {


    auto c     = read_access(state.c().asArray());
    auto del2f = write_access(laplacian_argument_->asArray());

    maDGForAll(i, 0, c.size(), { del2f[i] = pow(c[i], 3.0) - c[i]; });

    geom_.laplacian(*laplacian_argument_, *laplacian_rhs_term_);
  }

  /**
   * Compute biharmonic(c)
   */
  geom_.biharmonic(state.c(), *biharmonic_term_);

  /**
   * compute sigma*(c-m)
   */
  {
    auto          c     = read_access(state.c().asArray());
    auto          idx   = read_access(geom_.interiorIndices());
    auto          term  = write_access(linear_term_->asArray());
    const real_wp m     = this->m();
    const real_wp sigma = this->sigma();

    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      term[i]     = sigma * (c[i] - m);
    });
  }

  {
    ScalarSolutionState3D& dstate_dt = dynamic_cast<ScalarSolutionState3D&>(rhs);
    {
      auto r = write_access(dstate_dt.c().asArray());
      maDGForAll(i, 0, r.size(), { r[i] = real_wp(0); });
    }

    auto rhseval = write_access(dstate_dt.c().asArray());
    auto idx     = read_access(geom_.interiorIndices());
    auto del4    = read_access(biharmonic_term_->asArray());
    auto del2    = read_access(laplacian_rhs_term_->asArray());
    auto sigterm = read_access(linear_term_->asArray());

    const real_wp eps2 = this->eps2();

    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      rhseval[i]  = -eps2 * del4[i] + del2[i] - sigterm[i];
    });
  }
}

void
CahnHilliard3DFD::evalRHSFused(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const ScalarSolutionState3D& state     = dynamic_cast<const ScalarSolutionState3D&>(flovars);
  ScalarSolutionState3D&       dstate_dt = dynamic_cast<ScalarSolutionState3D&>(rhs);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));
  geom_.applyNeumannBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));

  /**
   * Compute  laplacian(u*c^3 -b*c)
   */
  {
    auto c     = read_access(state.c().asArray());
    auto del2f = write_access(laplacian_argument_->asArray());
    maDGForAll(i, 0, c.size(), { del2f[i] = pow(c[i], 3.0) - c[i]; });
  }

  {
    auto          laplacian_arg = read_access((*laplacian_argument_));
    auto          idx           = read_access(interiorIndices());
    const real_wp dx2           = geom_.dx() * geom_.dx();
    const real_wp dx4           = dx2 * dx2;
    auto          f             = read_access(state.c());

    auto rhseval = write_access(dstate_dt.c());

    const real_wp m     = this->m();
    const real_wp sigma = this->sigma();
    const real_wp eps2  = this->eps2();

    maDGForAll(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      laplacian_arg.getIJK(idx[ii], i, j, k);


      const real_wp del2 = (laplacian_arg(i - 1, j, k) + laplacian_arg(i + 1, j, k) + laplacian_arg(i, j - 1, k) +
                            laplacian_arg(i, j + 1, k) + laplacian_arg(i, j, k - 1) + laplacian_arg(i, j, k + 1) -
                            real_wp(6) * laplacian_arg(i, j, k)) /
                           dx2;

      const real_wp del4f =
          (real_wp(42) * f(i, j, k) -
           real_wp(12) *
               (f(i + 1, j, k) + f(i, j + 1, k) + f(i, j, k + 1) + f(i - 1, j, k) + f(i, j - 1, k) + f(i, j, k - 1)) +
           f(i + 2, j, k) + f(i, j + 2, k) + f(i, j, k + 2) + f(i - 2, j, k) + f(i, j - 2, k) + f(i, j, k - 2) +
           real_wp(2) * (f(i - 1, j - 1, k) + f(i - 1, j + 1, k) + f(i + 1, j - 1, k) + f(i + 1, j + 1, k) +
                         f(i - 1, j, k - 1) + f(i + 1, j, k - 1) + f(i, j - 1, k - 1) + f(i, j + 1, k - 1) +
                         f(i - 1, j, k + 1) + f(i + 1, j, k + 1) + f(i, j - 1, k + 1) + f(i, j + 1, k + 1))) /
          dx4;

      const real_wp linear_term = sigma * (f(i, j, k) - m);

      rhseval(i, j, k) = -eps2 * del4f + del2 - linear_term;
    });
  }
}


void
CahnHilliard3DFD::evalRHSFullyFused(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const ScalarSolutionState3D& state     = dynamic_cast<const ScalarSolutionState3D&>(flovars);
  ScalarSolutionState3D&       dstate_dt = dynamic_cast<ScalarSolutionState3D&>(rhs);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));
  geom_.applyNeumannBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));

  {
    auto          idx = read_access(interiorIndices());
    const real_wp dx2 = geom_.dx() * geom_.dx();
    const real_wp dx4 = dx2 * dx2;
    auto          f   = read_access(state.c());

    auto rhseval = write_access(dstate_dt.c());

    const real_wp m     = this->m();
    const real_wp sigma = this->sigma();
    const real_wp eps2  = this->eps2();

    maDGForAll(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      f.getIJK(idx[ii], i, j, k);

      const real_wp laparg_im1 = pow(f(i - 1, j, k), 3.0) - f(i - 1, j, k);
      const real_wp laparg_ip1 = pow(f(i + 1, j, k), 3.0) - f(i + 1, j, k);
      const real_wp laparg_jm1 = pow(f(i, j - 1, k), 3.0) - f(i, j - 1, k);
      const real_wp laparg_jp1 = pow(f(i, j + 1, k), 3.0) - f(i, j + 1, k);
      const real_wp laparg_km1 = pow(f(i, j, k - 1), 3.0) - f(i, j, k - 1);
      const real_wp laparg_kp1 = pow(f(i, j, k + 1), 3.0) - f(i, j, k + 1);
      const real_wp laparg_ijk = pow(f(i, j, k), 3.0) - f(i, j, k);

      const real_wp del2 =
          (laparg_im1 + laparg_ip1 + laparg_jm1 + laparg_jp1 + laparg_km1 + laparg_kp1 - real_wp(6) * laparg_ijk) / dx2;

      const real_wp del4f =
          (real_wp(42) * f(i, j, k) -
           real_wp(12) *
               (f(i + 1, j, k) + f(i, j + 1, k) + f(i, j, k + 1) + f(i - 1, j, k) + f(i, j - 1, k) + f(i, j, k - 1)) +
           f(i + 2, j, k) + f(i, j + 2, k) + f(i, j, k + 2) + f(i - 2, j, k) + f(i, j - 2, k) + f(i, j, k - 2) +
           real_wp(2) * (f(i - 1, j - 1, k) + f(i - 1, j + 1, k) + f(i + 1, j - 1, k) + f(i + 1, j + 1, k) +
                         f(i - 1, j, k - 1) + f(i + 1, j, k - 1) + f(i, j - 1, k - 1) + f(i, j + 1, k - 1) +
                         f(i - 1, j, k + 1) + f(i + 1, j, k + 1) + f(i, j - 1, k + 1) + f(i, j + 1, k + 1))) /
          dx4;

      const real_wp linear_term = sigma * (f(i, j, k) - m);

      rhseval(i, j, k) = -eps2 * del4f + del2 - linear_term;
    });
  }
}


std::unique_ptr<CSRMatrix>
CahnHilliard3DFD::createSparseMatrix() const
{
//  if (geom_.domain().xbc == BCType::Periodic && geom_.domain().ybc == BCType::Periodic &&
//      geom_.domain().zbc == BCType::Periodic)
//  {

    const int stencil_size = 25;
    const int npoints      = geom_.nInteriorPoints();
    const int nnz          = npoints * stencil_size;

    auto csr = std::make_unique<CSRMatrix>(npoints, npoints, nnz);

    {  // set all values to zero for sparsity pattern
      auto v = write_access(csr->values());
      maDGForAll(i, 0, v.size(), { v[i] = 0.0; });
    }

    {
      // same number of nonzeros per row so the row index is an equally strided array
      auto row_idx = write_access(csr->rowIdx());
      maDGForAll(irow, 0, row_idx.size(), { row_idx[irow] = stencil_size * irow; });
    }

    {  // set up column indexes
      const int ni = geom_.ni();
      const int nj = geom_.nj();
      const int nk = geom_.nk();

      auto periodic_idx = [&](int i, int j, int k) {
        if (i < 0) {
          i += ni;
        }
        if (j < 0) {
          j += nj;
        }
        if (k < 0) {
          k += nk;
        }
        if (i >= ni) {
          i -= ni;
        }
        if (j >= nj) {
          j -= nj;
        }
        if (k >= nk) {
          k -= nk;
        }
        return get1DindexFrom3D(i, j, k, 0, nj, nk);
      };

      int irow = 0;

      auto csr_row_idx = read_access_host(csr->rowIdx());
      auto csr_col_idx = write_access_host(csr->colIdx());

      for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
          for (int k = 0; k < nk; ++k) {
            std::vector<int_t> col_idx;
            col_idx.reserve(stencil_size);
            col_idx.push_back(periodic_idx(i, j, k));
            col_idx.push_back(periodic_idx(i + 1, j, k));
            col_idx.push_back(periodic_idx(i, j + 1, k));
            col_idx.push_back(periodic_idx(i, j, k + 1));
            col_idx.push_back(periodic_idx(i - 1, j, k));
            col_idx.push_back(periodic_idx(i, j - 1, k));
            col_idx.push_back(periodic_idx(i, j, k - 1));
            col_idx.push_back(periodic_idx(i + 2, j, k));
            col_idx.push_back(periodic_idx(i, j + 2, k));
            col_idx.push_back(periodic_idx(i, j, k + 2));
            col_idx.push_back(periodic_idx(i - 2, j, k));
            col_idx.push_back(periodic_idx(i, j - 2, k));
            col_idx.push_back(periodic_idx(i, j, k - 2));
            col_idx.push_back(periodic_idx(i - 1, j - 1, k));
            col_idx.push_back(periodic_idx(i - 1, j + 1, k));
            col_idx.push_back(periodic_idx(i + 1, j - 1, k));
            col_idx.push_back(periodic_idx(i + 1, j + 1, k));
            col_idx.push_back(periodic_idx(i - 1, j, k - 1));
            col_idx.push_back(periodic_idx(i + 1, j, k - 1));
            col_idx.push_back(periodic_idx(i, j - 1, k - 1));
            col_idx.push_back(periodic_idx(i, j + 1, k - 1));
            col_idx.push_back(periodic_idx(i - 1, j, k + 1));
            col_idx.push_back(periodic_idx(i + 1, j, k + 1));
            col_idx.push_back(periodic_idx(i, j - 1, k + 1));
            col_idx.push_back(periodic_idx(i, j + 1, k + 1));
            std::sort(col_idx.begin(), col_idx.end());
            col_idx.erase(std::unique(col_idx.begin(), col_idx.end()), col_idx.end());

            const int row_start = csr_row_idx[irow];
            const int row_end   = csr_row_idx[irow + 1];
            maDGForAllHost(icol, row_start, row_end, { csr_col_idx[icol] = col_idx[icol - row_start]; });
            irow++;
          }
        }
      }
    }

    return csr;
//  }
//  else {
//    Logger::get().FatalMessage("CahnHilliard3DFD::createSparseMatrix() only implemented for triply-periodic");
//  }

  Logger::get().FatalMessage("CahnHilliard3DFD::createSparseMatrix() not implemented");
  return nullptr;
}

namespace
{
MADG_HOST_DEVICE inline int periodic_idx(int i, int j, int k, int ni, int nj, int nk)
{
  if (i < 0) {
    i += ni;
  }
  if (j < 0) {
    j += nj;
  }
  if (k < 0) {
    k += nk;
  }
  if (i >= ni) {
    i -= ni;
  }
  if (j >= nj) {
    j -= nj;
  }
  if (k >= nk) {
    k -= nk;
  }
  return get1DindexFrom3D(i, j, k, 0, nj, nk);
}
}

std::vector<int> CahnHilliard3DFD::nNonZerosPerRow() const
{
//  Logger::get().FatalAssert(
//      geom_.domain().xbc == BCType::Periodic && geom_.domain().ybc == BCType::Periodic &&
//      geom_.domain().zbc == BCType::Periodic,
//      "CahnHilliard3DFD::evalJacobian() only supports triply periodic");

  const int stencil_size = 25;
  const int npoints      = geom_.nInteriorPoints();
  std::vector<int> nnz(npoints);
  std::fill(nnz.begin(), nnz.end(), stencil_size);
  return nnz;
}

void
CahnHilliard3DFD::evalJacobian(const SolutionState& flovars, double time, CSRMatrix& J)
{
  profile();
//  auto timer = Logger::get().timer("CahnHilliard3DFD::evalJacobian");

//  Logger::get().FatalAssert(
//      geom_.domain().xbc == BCType::Periodic && geom_.domain().ybc == BCType::Periodic &&
//          geom_.domain().zbc == BCType::Periodic,
//      "CahnHilliard3DFD::evalJacobian() only supports triply periodic");

  const auto& state = dynamic_cast<const ScalarSolutionState3D&>(flovars);

  auto f = read_access(state.c());

  const real_wp dx2   = geom_.dx() * geom_.dx();
  const real_wp dx4   = dx2 * dx2;
  const real_wp sigma = this->sigma();
  const real_wp eps2  = this->eps2();

  const int ni = geom_.ni();
  const int nj = geom_.nj();
  const int nk = geom_.nk();

//  auto periodic_idx = [&](int i, int j, int k) {
//    if (i < 0) {
//      i += ni;
//    }
//    if (j < 0) {
//      j += ni;
//    }
//    if (k < 0) {
//      k += ni;
//    }
//    if (i >= ni) {
//      i -= ni;
//    }
//    if (j >= ni) {
//      j -= ni;
//    }
//    if (k >= ni) {
//      k -= ni;
//    }
//    return get1DindexFrom3D(i, j, k, 0, ni, ni);
//  };


  auto col_idx = read_access(J.colIdx());
  auto row_idx = read_access(J.rowIdx());
  auto vals    = write_access(J.values());

  maDGForAll(i, 0, vals.size(), { vals[i] = 0; });

  maDGForAll(irow, 0, row_idx.size() - 1, {
    const int i = irow / (nj * nk);
    const int j = (irow / nk) % nj;
    const int k = irow % nk;

    const int row_start = row_idx[irow];
    const int row_end   = row_idx[irow + 1];

    for (int val_idx = row_start; val_idx < row_end; ++val_idx) {

      const int icol = col_idx[val_idx];

      /*
       * todo these ifs are stupid.
       * all the index offsets should be known in advance.
       */

      if (icol == irow)  // diagonal term
      {
        vals[val_idx] += -sigma;                                           // from -sigma*(c-m)
        vals[val_idx] += -eps2 * 42.0 / dx4;                               // from biharmonic
        vals[val_idx] += -6.0 * (3.0 * pow(f(i, j, k), 2.0) - 1.0) / dx2;  // from nonlinear laplacian
      }

      // local stencil terms covered by the laplacian
      if (icol == periodic_idx(i - 1, j, k, ni, nj, nk)) {
        vals[val_idx] += (3.0 * pow(f(i - 1, j, k), 2.0) - 1.0) / dx2;  // nonlinear laplacian
        vals[val_idx] += 12.0 * eps2 / dx4;                             // biharmonic
      }
      if (icol == periodic_idx(i + 1, j, k, ni, nj, nk)) {
        vals[val_idx] += (3.0 * pow(f(i + 1, j, k), 2.0) - 1.0) / dx2;  // nonlinear laplacian
        vals[val_idx] += 12.0 * eps2 / dx4;                             // biharmonic
      }
      if (icol == periodic_idx(i, j - 1, k, ni, nj, nk)) {
        vals[val_idx] += (3.0 * pow(f(i, j - 1, k), 2.0) - 1.0) / dx2;  // nonlinear laplacian
        vals[val_idx] += 12.0 * eps2 / dx4;                             // biharmonic
      }
      if (icol == periodic_idx(i, j + 1, k, ni, nj, nk)) {
        vals[val_idx] += (3.0 * pow(f(i, j + 1, k), 2.0) - 1.0) / dx2;  // nonlinear laplacian
        vals[val_idx] += 12.0 * eps2 / dx4;                             // biharmonic
      }
      if (icol == periodic_idx(i, j, k - 1, ni, nj, nk)) {
        vals[val_idx] += (3.0 * pow(f(i, j, k - 1), 2.0) - 1.0) / dx2;  // nonlinear laplacian
        vals[val_idx] += 12.0 * eps2 / dx4;                             // biharmonic
      }
      if (icol == periodic_idx(i, j, k + 1, ni, nj, nk)) {
        vals[val_idx] += (3.0 * pow(f(i, j, k + 1), 2.0) - 1.0) / dx2;  // nonlinear laplacian
        vals[val_idx] += 12.0 * eps2 / dx4;                             // biharmonic
      }

      // +/- 2 terms from biharmonic
      if (icol == periodic_idx(i + 2, j, k, ni, nj, nk) || icol == periodic_idx(i - 2, j, k, ni, nj, nk) || icol == periodic_idx(i, j + 2, k, ni, nj, nk) ||
          icol == periodic_idx(i, j - 2, k, ni, nj, nk) || icol == periodic_idx(i, j, k + 2, ni, nj, nk) || icol == periodic_idx(i, j, k - 2, ni, nj, nk)) {
        vals[val_idx] += -eps2 / dx4;
      }

      // diagonal neighbor terms from biharmonic
      if (icol == periodic_idx(i - 1, j - 1, k, ni, nj, nk) || icol == periodic_idx(i - 1, j + 1, k, ni, nj, nk) ||
          icol == periodic_idx(i + 1, j - 1, k, ni, nj, nk) || icol == periodic_idx(i + 1, j + 1, k, ni, nj, nk) ||

          icol == periodic_idx(i - 1, j, k - 1, ni, nj, nk) || icol == periodic_idx(i + 1, j, k - 1, ni, nj, nk) ||
          icol == periodic_idx(i, j - 1, k - 1, ni, nj, nk) || icol == periodic_idx(i, j + 1, k - 1, ni, nj, nk) ||

          icol == periodic_idx(i - 1, j, k + 1, ni, nj, nk) || icol == periodic_idx(i + 1, j, k + 1, ni, nj, nk) ||
          icol == periodic_idx(i, j - 1, k + 1, ni, nj, nk) || icol == periodic_idx(i, j + 1, k + 1, ni, nj, nk)) {
        vals[val_idx] += -eps2 * 2.0 / dx4;
      }
    }
  });

//  Logger::get().TraceMessage("Diagonal dominance without time contribution: " + std::to_string(J.diagonalDominance()));
}


std::unique_ptr<TimeIntegrableRHS> CahnHilliard3DFD::clone(SpatialDiscretization& geom) const
{
  return std::make_unique<CahnHilliard3DFD>(dynamic_cast<Discretization3DCart&>(geom), this->params());
}
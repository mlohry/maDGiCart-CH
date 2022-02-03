#pragma once

#include "spatial_discretization/spatial_discretization.hpp"

#include <math.h>

inline CartesianDomainDefinition
generateTriplyPeriodicDomain(int nx)
{
  CartesianDomainDefinition domain;
  domain.nx    = nx;
  domain.ny    = nx;
  domain.nz    = nx;
  domain.xbeg  = -M_PI;
  domain.xend  = M_PI;
  domain.ybeg  = -M_PI;
  domain.zbeg  = -M_PI;
  domain.xbc   = BCType::Periodic;
  domain.ybc   = BCType::Periodic;
  domain.zbc   = BCType::Periodic;
  domain.nhalo = 2;
  return domain;
}

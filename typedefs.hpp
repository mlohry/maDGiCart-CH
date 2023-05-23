#pragma once

#include <cstdint>


using int_t  = std::int32_t;

#ifdef MADG_USE_SINGLE_PRECISION
using real_t = float;
using real_wp = float;
#else
using real_t = double;
using real_wp = double;
#endif

#include <gtest/gtest.h>

#include "data_structures/exec_includes.hpp"
#include "data_structures/managed_array.hpp"

#include <iostream>


TEST(ManagedArray, BoundsCheckAssertionDebug)
{
#ifndef NDEBUG
  const int             N = 10;
  ManagedArray<real_wp> array({}, "TestArray", N);
  auto                  array_read = read_access_host(array);
  ASSERT_DEATH(
      {
        const real_wp x = array_read[-1];
      },
      "");
  ASSERT_DEATH(
      {
        const real_wp x = array_read[N];
      },
      "");

#ifdef __NVCC__
  auto array_read_device = read_access(array);
  ASSERT_DEATH(
      {
        const real_wp x = array_read_device[-1];
      },
      "");
  ASSERT_DEATH(
      {
        const real_wp x = array_read_device[N];
      },
      "");
#endif

#else
  std::cout << "Unit test " << ::testing::UnitTest::GetInstance()->current_test_info()->name()
            << " only valid in debug mode.\n";
#endif
}

#include <gtest/gtest.h>

#include "typedefs.hpp"

#include "data_structures/exec_includes.hpp"
#include "testing/utils/order_of_accuracy.hpp"
#include "data_structures/managed_array_3d.hpp"


TEST(Discretiztion3DCart, DataStructures)
{
  const int N = 2;
  const int Nhalo = 1;
  ManagedArray3DOwning<double> array3d({}, "arrayname", N, N, N, Nhalo);

  auto arr = read_access_host(array3d);

  int idx1d = 0;
  for (int i = -Nhalo; i < N+Nhalo; ++i){
    for (int j = -Nhalo; j < N+Nhalo; ++j){
      for (int k = -Nhalo; k < N+Nhalo; ++k){
        std::cout << "i,j,k: " << i << " " << j << " " << k << " idx1d: " << idx1d << "\n";

        int ii, jj, kk;
        arr.getIJK(idx1d, ii, jj, kk);
//        std::cout << "ii,jj,kk: " << ii << " " << jj << " " << kk << " idx1d: " << idx1d << "\n";
        ASSERT_EQ(i, ii);
        ASSERT_EQ(j, jj);
        ASSERT_EQ(k, kk);

        idx1d++;
      }
    }
  }

}

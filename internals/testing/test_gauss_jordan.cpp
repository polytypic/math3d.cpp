#include "math3d_v1/mtx.hpp"
#include "math3d_v1/transform.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;

auto test_gauss_jordan = test([]() {
  auto a_prime = mtx<float, 3>{{{3, 1, 4}, {1, 5, 9}, {2, 6, 5}}};
  auto a = a_prime;
  auto y = make_identity<float, 3>();

  gauss_jordan(&a, &y);

  auto identity = y * a_prime;

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      verify(std::abs(identity[i][j] - (i == j)) < 0.01);
});

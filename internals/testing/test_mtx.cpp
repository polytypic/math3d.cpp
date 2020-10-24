#include "math3d_v1/mtx.hpp"
#include "math3d_v1/transform.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;

template <class S, size_t N>
static bool is_identity(const mtx<S, N> &m, S epsilon = 0.01f) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (epsilon < std::abs(m[i][j] - (i == j)))
        return false;
  return true;
}

auto test_gauss_jordan = test([]() {
  auto a_prime = mtx<float, 3>{{{3, 1, 4}, {1, 5, 9}, {2, 6, 5}}};
  auto a = a_prime;
  auto y = make_identity<float, 3>();

  gauss_jordan(&a, &y);

  verify(is_identity(y * a_prime));
});

auto test_inverse = test([]() {
  auto a = mtx<float, 3>{{{1, 4, 1}, {5, 9, 2}, {6, 5, 3}}};

  verify(is_identity(a * inverse(a)));
});

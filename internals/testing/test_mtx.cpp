#include "math3d_v1/mtx.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;

template <class S, size_t N>
static bool is_identity(const mtx<S, N> &m, S epsilon = 0.01f) {
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j)
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

auto test_ops = test([]() {
  const mtx<float, 1> m3{{{3.0f}}};
  const mtx<float, 1> m6{{{6.0f}}};
  verify((m3 + m6)[0][0] == 9.0f);
  verify((m3 - m6)[0][0] == -3.0f);
});

auto test_transpose = test([]() {
  const auto m = from_columns(make_vec(1, 2, 3), make_vec(4, 5, 6));
  const auto t = transpose(m);
  verify(m.columns() == t.rows());
  verify(m.rows() == t.columns());
  for (size_t i = 0; i < m.rows(); ++i)
    for (size_t j = 0; j < m.columns(); ++j)
      verify(m[i][j] == t[j][i]);
});

auto test_make_identity =
    test([]() { verify(is_identity(make_identity<float, 4>())); });

auto test_from_diagonal = test([]() {
  verify(is_identity(from_diagonal(make_vec(2.0f, 2.0f)) *
                     from_diagonal(make_vec(0.5f, 0.5f))));
});

auto test_mtx_vec = test([]() {
  auto v = from_diagonal(make_vec(2.0f, 4.0f)) * make_vec(0.5f, 0.25f);
  verify(v[0] == 1 && v[1] == 1);
});

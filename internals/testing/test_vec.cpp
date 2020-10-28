#include "math3d_v1/vec.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;

auto test_mag = test([]() { verify(mag(vec(vec(1, 0), 0)) == 1); });

auto test_lerp = test([]() { verify(lerp(vec(1), vec(0), 0.5f)[0] == 0.5f); });

auto test_normalize = test([]() { verify(mag(normalize(vec(0, 2, 0))) == 1); });

auto test_zero_vec = test([]() { verify(norm(zero_vec<float, 4>()) == 0); });

auto test_ops = test([]() {
  verify((vec(1) + vec(2))[0] == 3);
  verify((vec(1) - vec(2))[0] == -1);
  verify((vec(2) * vec(3))[0] == 6);
  verify((vec(6) / vec(3))[0] == 2);
  verify((vec(3) % vec(2))[0] == 1);

  verify((vec(1) + 2)[0] == 3);
  verify((vec(1) - 2)[0] == -1);
  verify((vec(2) * 3)[0] == 6);
  verify((vec(6) / 3)[0] == 2);
  verify((vec(3) % 2)[0] == 1);

  verify((1 + vec(2))[0] == 3);
  verify((1 - vec(2))[0] == -1);
  verify((2 * vec(3))[0] == 6);
  verify((6 / vec(3))[0] == 2);
  verify((3 % vec(2))[0] == 1);
});

auto test_cross =
    test([]() { verify(cross(vec(1, 0, 0), vec(0, 1, 0))[2] == 1); });

auto test_homogenize =
    test([]() { verify(homogenize(vec(6.0f, 3.0f))[0] == 2.0f); });

auto test_sub =
    test([]() { verify((sub<2>(vec(vec(1), 2, 3))).dimensions() == 2); });

#include "math3d_v1/vec.hpp"

#include "data_v1/array.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;
using namespace data_v1;

auto test_mag = test([]() { verify(mag(make_vec(1, 0, 0)) == 1); });

auto test_lerp =
    test([]() { verify(lerp(make_vec(1), make_vec(0), 0.5f)[0] == 0.5f); });

auto test_normalize =
    test([]() { verify(mag(normalize(make_vec(0, 2, 0))) == 1); });

auto test_zero_vec = test([]() { verify(norm(zero_vec<float, 4>()) == 0); });

auto test_ops = test([]() {
  verify((make_vec(1) + make_vec(2))[0] == 3);
  verify((make_vec(1) - make_vec(2))[0] == -1);
  verify((make_vec(2) * make_vec(3))[0] == 6);
  verify((make_vec(6) / make_vec(3))[0] == 2);
  verify((make_vec(3) % make_vec(2))[0] == 1);

  verify((make_vec(1) + 2)[0] == 3);
  verify((make_vec(1) - 2)[0] == -1);
  verify((make_vec(2) * 3)[0] == 6);
  verify((make_vec(6) / 3)[0] == 2);
  verify((make_vec(3) % 2)[0] == 1);

  verify((1 + make_vec(2))[0] == 3);
  verify((1 - make_vec(2))[0] == -1);
  verify((2 * make_vec(3))[0] == 6);
  verify((6 / make_vec(3))[0] == 2);
  verify((3 % make_vec(2))[0] == 1);
});

auto test_cross =
    test([]() { verify(cross(make_vec(1, 0, 0), make_vec(0, 1, 0))[2] == 1); });

auto test_homogenize =
    test([]() { verify(homogenize(make_vec(6.0f, 3.0f))[0] == 2.0f); });

auto test_sub =
    test([]() { verify(size(sub<2>(make_vec(1, 2, 3)).values) == 2); });

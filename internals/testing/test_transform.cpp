#include "math3d_v1/transform.hpp"
#include "math3d_v1/trig.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;

auto test_scaling = test([]() {
  auto v = scaling(vec(2.0f, 4.0f, 1.0f)) * vec(0.5f, 0.25f, 1.0f, 0.0f);
  verify(v[0] == 1 && v[1] == 1 && v[2] == 1);
});

auto test_translation = test([]() {
  auto v = translation(vec(2.0f, 4.0f, 1.0f)) * vec(-1.0f, -3.0f, 0.0f, 1.0f);
  verify(v[0] == 1 && v[1] == 1 && v[2] == 1);
});

auto test_rotation = test([]() {
  for (size_t i = 0; i < 3; ++i) {
    auto m = rotation<float>(float(pi) * 0.5f, i);
    auto u = zero_vec<float, 4>();
    u[(i + 1) % 3] = 1;
    auto v = m * u;
    verify(std::abs(v[(i + 0) % 3] - 0) < 0.001f);
    verify(std::abs(v[(i + 1) % 3] - 0) < 0.001f);
    verify(std::abs(v[(i + 2) % 3] - 1) < 0.001f);
  }
});

auto test_perspective = test([]() {
  auto fov = from_angle(90.0f);

  auto near = 1.0f;
  auto far = 5.0f;

  auto m = perspective(fov, near, far);

  {
    auto u = vec(near * std::tan(fov * 0.5f), 0.0f, -near, 1.0f);

    auto v = homogenize(m * u);

    verify(std::abs(v[2]) < 0.001f);
    verify(std::abs(v[0] - 1) < 0.001f);
  }

  {
    auto u = vec(0.0f, far * std::tan(fov * 0.5f), -far, 1.0f);

    auto v = homogenize(m * u);

    verify(std::abs(v[2] - 1) < 0.001f);
    verify(std::abs(v[1] - 1) < 0.001f);
  }
});

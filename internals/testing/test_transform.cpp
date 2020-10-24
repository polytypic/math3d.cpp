#include "math3d_v1/transform.hpp"

#include "testing_v1/test.hpp"

using namespace testing_v1;
using namespace math3d_v1;

auto test_mtx_vec = test([]() {
  auto v = make_scaling(make_vec(2.0f, 4.0f, 1.0f)) *
           make_vec(0.5f, 0.25f, 1.0f, 0.0f);
  verify(v[0] == 1 && v[1] == 1 && v[2] == 1);
});

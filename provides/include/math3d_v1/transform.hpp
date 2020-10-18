#pragma once

#include "math3d_v1/synopsis.hpp"

#include <cassert>

template <class Scalar>
auto math3d_v1::make_projection(Scalar fov, Scalar near, Scalar far) {
  static const Scalar fov_scale = 0.5 * pi / 180;
  Scalar scale = 1 / tan(fov * fov_scale);
  Scalar near_far_sub_inv = 1 / (near - far);
  return mtx<Scalar, 4>{
      {{scale, 0, 0, 0},
       {0, scale, 0, 0},
       {0, 0, far * near_far_sub_inv, far * near_far_sub_inv * near},
       {0, 0, -1, 0}}};
}

template <class Scalar>
auto math3d_v1::make_translation(const vec<Scalar, 3> &trans) {
  return mtx<Scalar, 4>{{{1, 0, 0, trans[0]},
                         {0, 1, 0, trans[1]},
                         {0, 0, 1, trans[2]},
                         {0, 0, 0, 1}}};
}

template <class Scalar>
auto math3d_v1::make_scaling(const vec<Scalar, 3> &scale) {
  return mtx<Scalar, 4>{{{scale[0], 0, 0, 0},
                         {0, scale[1], 0, 0},
                         {0, 0, scale[2], 0},
                         {0, 0, 0, 1}}};
}

template <class Scalar>
auto math3d_v1::make_rotation(Scalar radians, size_t axis) {
  assert(axis < 3);

  mtx<Scalar, 4> result{};

  result[axis][axis] = 1;
  result[3][3] = 1;

  size_t axis1 = (axis + 1) % 3;
  size_t axis2 = (axis + 2) % 3;

  Scalar c = cos(radians);
  Scalar s = sin(radians);

  result[axis1][axis1] = c;
  result[axis1][axis2] = -s;

  result[axis2][axis1] = s;
  result[axis2][axis2] = c;

  return result;
}
#pragma once

#include "math3d_v1/synopsis.hpp"

template <class Scalar> inline auto math3d_v1::from_angle(Scalar angle) {
  return angle * Scalar(pi / 180);
}

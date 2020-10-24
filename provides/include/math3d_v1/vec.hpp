#pragma once

#include "math3d_v1/synopsis.hpp"

#include "math3d_v1/ops.hpp"

#include <cassert>
#include <cmath>

template <class S, size_t N> auto math3d_v1::vec<S, N>::dimensions() const {
  return N;
}

template <class S, size_t N>
const auto &math3d_v1::vec<S, N>::operator[](size_t i) const {
  assert(i < N);
  return values[i];
}

template <class S, size_t N> auto &math3d_v1::vec<S, N>::operator[](size_t i) {
  assert(i < N);
  return values[i];
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::bin_op(const vec<SL, N> &lhs, const vec<SR, N> &rhs) {
  vec<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs[i], rhs[i]);
  return result;
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::bin_op(SL lhs, const vec<SR, N> &rhs) {
  vec<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs, rhs[i]);
  return result;
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::bin_op(const vec<SL, N> &lhs, SR rhs) {
  vec<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs[i], rhs);
  return result;
}

#define MAKE(name, op)                                                         \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(const vec<SL, N> &lhs, const vec<SR, N> &rhs) {           \
    return bin_op<name>(lhs, rhs);                                             \
  }                                                                            \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(SL lhs, const vec<SR, N> &rhs) {                          \
    return bin_op<name>(lhs, rhs);                                             \
  }                                                                            \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(const vec<SL, N> &lhs, SR rhs) {                          \
    return bin_op<name>(lhs, rhs);                                             \
  }

MAKE(add_op, operator+)
MAKE(sub_op, operator-)
MAKE(mul_op, operator*)
MAKE(div_op, operator/)
MAKE(rem_op, operator%)
#undef MAKE

template <class SL, class SR, size_t N>
auto math3d_v1::dot(const vec<SL, N> &lhs, const vec<SR, N> &rhs) {
  auto result = lhs[0] * rhs[0];
  for (size_t i = 1; i < N; ++i)
    result += lhs[i] * rhs[i];
  return result;
}

template <class S, size_t N> auto math3d_v1::norm(const vec<S, N> &v) {
  return dot(v, v);
}

template <class S, size_t N> auto math3d_v1::mag(const vec<S, N> &v) {
  return std::sqrt(norm(v));
}

template <class S, size_t N> auto math3d_v1::zero_vec() { return vec<S, N>{}; }

template <class SL, class SR, class ST, size_t N>
auto math3d_v1::lerp(const vec<SL, N> &lhs, const vec<SR, N> &rhs, ST t) {
  vec<decltype(SL() * ST() + SR() * ST()), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = t * (rhs[i] - lhs[i]) + lhs[i];
  return result;
}

template <class SL, class SR>
auto math3d_v1::cross(const vec<SL, 3> &lhs, const vec<SR, 3> &rhs) {
  return make_vec(lhs[1] * rhs[2] - lhs[2] * rhs[1],
                  lhs[2] * rhs[0] - lhs[0] * rhs[2],
                  lhs[0] * rhs[1] - lhs[1] * rhs[0]);
}

template <class S, size_t N> auto math3d_v1::normalize(const vec<S, N> &v) {
  return v * (1 / mag(v));
}

template <class S, size_t N> auto math3d_v1::homogenize(const vec<S, N> &v) {
  return v * (1 / v[N - 1]);
}

template <size_t N, size_t I, class S, size_t M>
auto math3d_v1::sub(const vec<S, M> &v) {
  static_assert(I + N <= M);
  vec<S, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = v[i + I];
  return result;
}

template <class Scalar, class... Scalars>
auto math3d_v1::make_vec(Scalar value, Scalars... values) {
  return vec<decltype(+(value + ... + values)), 1 + sizeof...(Scalars)>{
      value, values...};
}

template <class Scalar, size_t N, class... Scalars>
auto math3d_v1::make_vec(const vec<Scalar, N> &v, Scalars... values) {
  vec<decltype(+(Scalar() + ... + values)), N + sizeof...(Scalars)> result;
  size_t i = 0;
  do {
    result[i] = v[i];
  } while (++i < N);
  ((result[i++] = values), ...);
  return result;
}

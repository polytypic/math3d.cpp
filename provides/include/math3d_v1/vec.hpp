#pragma once

#include "math3d_v1/synopsis.hpp"

#include "math3d_v1/ops.hpp"

#include <cassert>
#include <cmath>

template <class S, size_t N> auto math3d_v1::vec_t<S, N>::dimensions() const {
  return N;
}

template <class S, size_t N>
const auto &math3d_v1::vec_t<S, N>::operator[](size_t i) const {
  assert(i < N);
  return values[i];
}

template <class S, size_t N>
auto &math3d_v1::vec_t<S, N>::operator[](size_t i) {
  assert(i < N);
  return values[i];
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::bin_op(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs) {
  vec_t<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs[i], rhs[i]);
  return result;
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::bin_op(SL lhs, const vec_t<SR, N> &rhs) {
  vec_t<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs, rhs[i]);
  return result;
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::bin_op(const vec_t<SL, N> &lhs, SR rhs) {
  vec_t<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs[i], rhs);
  return result;
}

#define MAKE(name, op)                                                         \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs) {       \
    return bin_op<name>(lhs, rhs);                                             \
  }                                                                            \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(SL lhs, const vec_t<SR, N> &rhs) {                        \
    return bin_op<name>(lhs, rhs);                                             \
  }                                                                            \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(const vec_t<SL, N> &lhs, SR rhs) {                        \
    return bin_op<name>(lhs, rhs);                                             \
  }

MAKE(add_op_t, operator+)
MAKE(sub_op_t, operator-)
MAKE(mul_op_t, operator*)
MAKE(div_op_t, operator/)
MAKE(rem_op_t, operator%)
#undef MAKE

template <class SL, class SR, size_t N>
auto math3d_v1::dot(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs) {
  auto result = lhs[0] * rhs[0];
  for (size_t i = 1; i < N; ++i)
    result += lhs[i] * rhs[i];
  return result;
}

template <class S, size_t N> auto math3d_v1::norm(const vec_t<S, N> &v) {
  return dot(v, v);
}

template <class S, size_t N> auto math3d_v1::mag(const vec_t<S, N> &v) {
  return std::sqrt(norm(v));
}

template <class S, size_t N> auto math3d_v1::zero_vec() {
  return vec_t<S, N>{};
}

template <class SL, class SR, class ST, size_t N>
auto math3d_v1::lerp(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs, ST t) {
  vec_t<decltype(SL() * ST() + SR() * ST()), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = t * (rhs[i] - lhs[i]) + lhs[i];
  return result;
}

template <class SL, class SR>
auto math3d_v1::cross(const vec_t<SL, 3> &lhs, const vec_t<SR, 3> &rhs) {
  return vec(lhs[1] * rhs[2] - lhs[2] * rhs[1],
             lhs[2] * rhs[0] - lhs[0] * rhs[2],
             lhs[0] * rhs[1] - lhs[1] * rhs[0]);
}

template <class S, size_t N> auto math3d_v1::normalize(const vec_t<S, N> &v) {
  return v * (1 / mag(v));
}

template <class S, size_t N> auto math3d_v1::homogenize(const vec_t<S, N> &v) {
  return v * (1 / v[N - 1]);
}

template <size_t N, size_t I, class S, size_t M>
auto math3d_v1::sub(const vec_t<S, M> &v) {
  static_assert(I + N <= M);
  vec_t<S, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = v[i + I];
  return result;
}

template <class Scalar, class... Scalars>
auto math3d_v1::vec(Scalar value, Scalars... values) {
  return vec_t<decltype(+(value + ... + values)), 1 + sizeof...(Scalars)>{
      value, values...};
}

template <class Scalar, size_t N, class... Scalars>
auto math3d_v1::vec(const vec_t<Scalar, N> &v, Scalars... values) {
  vec_t<decltype(+(Scalar() + ... + values)), N + sizeof...(Scalars)> result;
  size_t i = 0;
  do {
    result[i] = v[i];
  } while (++i < N);
  ((result[i++] = values), ...);
  return result;
}

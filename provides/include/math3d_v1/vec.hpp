#pragma once

#include "math3d_v1/synopsis.hpp"

#include "math3d_v1/ops.hpp"

#include <cassert>

template <class S, size_t N>
const S &math3d_v1::vec<S, N>::operator[](size_t i) const {
  assert(i < N);
  return values[i];
}

template <class S, size_t N> S &math3d_v1::vec<S, N>::operator[](size_t i) {
  assert(i < N);
  return values[i];
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::binOp(const vec<SL, N> &lhs, const vec<SR, N> &rhs) {
  vec<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs[i], rhs[i]);
  return result;
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::binOp(SL lhs, const vec<SR, N> &rhs) {
  vec<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs, rhs[i]);
  return result;
}

template <class BinOp, class SL, class SR, size_t N>
auto math3d_v1::binOp(const vec<SL, N> &lhs, SR rhs) {
  vec<decltype(BinOp::apply(SL(), SR())), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = BinOp::apply(lhs[i], rhs);
  return result;
}

#define MAKE(name, op)                                                         \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(const vec<SL, N> &lhs, const vec<SR, N> &rhs) {           \
    return binOp<name>(lhs, rhs);                                              \
  }                                                                            \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(SL lhs, const vec<SR, N> &rhs) {                          \
    return binOp<name>(lhs, rhs);                                              \
  }                                                                            \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto math3d_v1::op(const vec<SL, N> &lhs, SR rhs) {                          \
    return binOp<name>(lhs, rhs);                                              \
  }

MAKE(add_op, operator+)
MAKE(sub_op, operator-)
MAKE(mul_op, operator*)
MAKE(div_op, operator/)
MAKE(rem_op, operator%)
#undef MAKE

template <class SL, class SR, size_t N>
auto math3d_v1::dot(const vec<SL, N> &lhs, const vec<SR, N> &rhs) {
  decltype(SL() * SR()) result = lhs[0] * rhs[0];
  for (size_t i = 1; i < N; ++i)
    result += lhs[i] * rhs[i];
  return result;
}

template <class S, size_t N> auto math3d_v1::norm(const vec<S, N> &v) {
  return dot(v, v);
}

template <class S, size_t N> auto math3d_v1::mag(const vec<S, N> &v) {
  return sqrt(norm(v));
}

template <class S, size_t N> math3d_v1::vec<S, N> math3d_v1::zero_vec() {
  return {};
}

template <class SL, class SR, class ST, size_t N>
auto math3d_v1::lerp(const vec<SL, N> &lhs, const vec<SR, N> &rhs, ST t) {
  vec<decltype(SL() * ST() + SR() * ST()), N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = t * (rhs[i] - lhs[i]) + lhs[i];
  return result;
}
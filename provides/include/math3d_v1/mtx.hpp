#pragma once

#include "math3d_v1/ops.hpp"
#include "math3d_v1/synopsis.hpp"

#include <cassert>

template <class Scalar, size_t R, size_t C>
const Scalar (&math3d_v1::mtx<Scalar, R, C>::operator[](size_t i) const)[C] {
  assert(i < R);
  return values[i];
}

template <class Scalar, size_t R, size_t C>
Scalar (&math3d_v1::mtx<Scalar, R, C>::operator[](size_t i))[C] {
  assert(i < R);
  return values[i];
}

template <class BinOp, class SL, class SR, size_t R, size_t C>
auto math3d_v1::binOp(const mtx<SL, R, C> &lhs, const mtx<SR, R, C> &rhs) {
  mtx<decltype(BinOp::apply(SL(), SR())), R, C> result;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j)
      result[i][j] = BinOp::apply(lhs[i][j], rhs[i][j]);
  return result;
}

#define BIN_OP(name, op)                                                       \
  template <class SL, class SR, size_t R, size_t C>                            \
  auto math3d_v1::op(const mtx<SL, R, C> &lhs, const mtx<SR, R, C> &rhs) {     \
    return binOp<name>(lhs, rhs);                                              \
  }

BIN_OP(add_op, operator+)
BIN_OP(sub_op, operator-)
#undef BIN_OP

template <class SL, class SR, size_t R, size_t RC, size_t C>
auto math3d_v1::operator*(const mtx<SL, R, RC> &lhs,
                          const mtx<SR, RC, C> &rhs) {
  using S = decltype(SL() * SR());
  mtx<S, R, C> result;
  for (size_t i = 0; i < R; ++i) {
    for (size_t j = 0; j < C; ++j) {
      S accum = lhs[i][0] * rhs[0][j];
      for (size_t k = 1; k < RC; ++k)
        accum += lhs[i][k] * rhs[k][j];
      result[i][j] = accum;
    }
  }
  return result;
}

template <class S, size_t R, size_t C>
auto math3d_v1::transpose(const mtx<S, R, C> &m) {
  mtx<S, C, R> result;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j)
      result[j][i] = m[i][j];
  return result;
}

template <class S, size_t N> auto math3d_v1::set_identity(mtx<S, N> *out) {
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j)
      (*out)[i][j] = +(i == j);
}

template <class S, size_t N> auto math3d_v1::make_identity() {
  mtx<S, N> result;
  set_identity(&result);
  return result;
}

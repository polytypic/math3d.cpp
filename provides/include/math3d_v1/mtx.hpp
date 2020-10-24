#pragma once

#include "math3d_v1/ops.hpp"
#include "math3d_v1/synopsis.hpp"

#include <cassert>
#include <cmath>

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
auto math3d_v1::bin_op(const mtx<SL, R, C> &lhs, const mtx<SR, R, C> &rhs) {
  mtx<decltype(BinOp::apply(SL(), SR())), R, C> result;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j)
      result[i][j] = BinOp::apply(lhs[i][j], rhs[i][j]);
  return result;
}

#define BIN_OP(name, op)                                                       \
  template <class SL, class SR, size_t R, size_t C>                            \
  auto math3d_v1::op(const mtx<SL, R, C> &lhs, const mtx<SR, R, C> &rhs) {     \
    return bin_op<name>(lhs, rhs);                                             \
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
  assert(out);
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j)
      (*out)[i][j] = S(i == j);
}

template <class S, size_t N> auto math3d_v1::make_identity() {
  mtx<S, N> result;
  set_identity(&result);
  return result;
}

template <class S, size_t N, size_t M>
auto math3d_v1::gauss_jordan(mtx<S, N> *pa, mtx<S, N, M> *py) {
  assert(pa && py);
  auto &a = *pa;
  auto &y = *py;
  for (size_t column = 0; column < N; ++column) {
    {
      size_t max_row = column;
      for (size_t row = column + 1; row < N; ++row) {
        if (std::abs(a[row][column]) > std::abs(a[max_row][column]))
          max_row = row;
      }

      for (size_t j = 0; j < N; ++j)
        std::swap(a[max_row][j], a[column][j]);
      for (size_t j = 0; j < M; ++j)
        std::swap(y[max_row][j], y[column][j]);
    }

    {
      auto m = 1 / a[column][column];
      for (size_t j = column; j < N; ++j)
        a[column][j] *= m;
      for (size_t j = 0; j < M; ++j)
        y[column][j] *= m;
    }

    for (size_t i = 0; i < N; ++i) {
      if (i != column) {
        auto m = a[i][column];
        for (size_t j = column; j < N; ++j)
          a[i][j] -= a[column][j] * m;
        for (size_t j = 0; j < M; ++j)
          y[i][j] -= y[column][j] * m;
      }
    }
  }
}

template <class S, size_t N> auto math3d_v1::inverse(const mtx<S, N> &m) {
  auto system = m;
  auto result = make_identity<S, N>();
  gauss_jordan(&system, &result);
  return result;
}

template <class SL, size_t N, size_t M, class SR>
auto math3d_v1::operator*(const mtx<SL, N, M> &m, const vec<SR, M> &v) {
  using S = decltype(SL() * SR());
  vec<S, N> result;
  for (size_t i = 0; i < N; ++i) {
    S accum = m[i][0] * v[0];
    for (size_t j = 1; j < M; ++j) {
      accum += m[i][j] * v[j];
    }
    result[i] = accum;
  }
  return result;
}

template <size_t Rows, class... Scalars>
auto math3d_v1::from_columns(const vec<Scalars, Rows> &... columns) {
  using Scalar = decltype((Scalars() + ...));
  mtx<Scalar, Rows, sizeof...(Scalars)> result;
  for (size_t i = 0; i < Rows; ++i) {
    size_t j = 0;
    ((result[i][j++] = columns[i]), ...);
  }
  return result;
}

template <class Scalar, size_t Rows>
auto math3d_v1::from_diagonal(const vec<Scalar, Rows> &diagonal) {
  mtx<Scalar, Rows> result{};
  for (size_t i = 0; i < Rows; ++i)
    result[i][i] = diagonal[i];
  return result;
}

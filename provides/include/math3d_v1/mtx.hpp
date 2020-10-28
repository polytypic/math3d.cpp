#pragma once

#include "math3d_v1/vec.hpp"

template <class Scalar, size_t R, size_t C>
auto math3d_v1::mtx_t<Scalar, R, C>::rows() const {
  return R;
}

template <class Scalar, size_t R, size_t C>
auto math3d_v1::mtx_t<Scalar, R, C>::columns() const {
  return C;
}

template <class Scalar, size_t R, size_t C>
const auto &math3d_v1::mtx_t<Scalar, R, C>::operator[](size_t i) const {
  assert(i < R);
  return values[i];
}

template <class Scalar, size_t R, size_t C>
auto &math3d_v1::mtx_t<Scalar, R, C>::operator[](size_t i) {
  assert(i < R);
  return values[i];
}

template <class BinOp, class SL, class SR, size_t R, size_t C>
auto math3d_v1::bin_op(const mtx_t<SL, R, C> &lhs, const mtx_t<SR, R, C> &rhs) {
  mtx_t<decltype(BinOp::apply(SL(), SR())), R, C> result;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j)
      result[i][j] = BinOp::apply(lhs[i][j], rhs[i][j]);
  return result;
}

#define MAKE(name, op)                                                         \
  template <class SL, class SR, size_t R, size_t C>                            \
  auto math3d_v1::op(const mtx_t<SL, R, C> &lhs, const mtx_t<SR, R, C> &rhs) { \
    return bin_op<name>(lhs, rhs);                                             \
  }

MAKE(add_op_t, operator+)
MAKE(sub_op_t, operator-)
#undef MAKE

template <class SL, class SR, size_t R, size_t RC, size_t C>
auto math3d_v1::operator*(const mtx_t<SL, R, RC> &lhs,
                          const mtx_t<SR, RC, C> &rhs) {
  using S = decltype(SL() * SR());
  mtx_t<S, R, C> result;
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
auto math3d_v1::transpose(const mtx_t<S, R, C> &m) {
  mtx_t<S, C, R> result;
  for (size_t i = 0; i < R; ++i)
    for (size_t j = 0; j < C; ++j)
      result[j][i] = m[i][j];
  return result;
}

template <class S, size_t N> auto math3d_v1::set_identity(mtx_t<S, N> *out) {
  assert(out);
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j)
      (*out)[i][j] = S(i == j);
}

template <class S, size_t N> auto math3d_v1::identity() {
  mtx_t<S, N> result;
  set_identity(&result);
  return result;
}

template <class S, size_t N, size_t M>
auto math3d_v1::gauss_jordan(mtx_t<S, N> *pa, mtx_t<S, N, M> *py) {
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

template <class S, size_t N> auto math3d_v1::inverse(const mtx_t<S, N> &m) {
  auto system = m;
  auto result = identity<S, N>();
  gauss_jordan(&system, &result);
  return result;
}

template <class SL, size_t N, size_t M, class SR>
auto math3d_v1::operator*(const mtx_t<SL, N, M> &m, const vec_t<SR, M> &v) {
  using S = decltype(SL() * SR());
  vec_t<S, N> result;
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
auto math3d_v1::from_columns(const vec_t<Scalars, Rows> &...columns) {
  using Scalar = decltype((Scalars() + ...));
  mtx_t<Scalar, Rows, sizeof...(Scalars)> result;
  for (size_t i = 0; i < Rows; ++i) {
    size_t j = 0;
    ((result[i][j++] = columns[i]), ...);
  }
  return result;
}

template <class Scalar, size_t Rows>
auto math3d_v1::from_diagonal(const vec_t<Scalar, Rows> &diagonal) {
  mtx_t<Scalar, Rows> result{};
  for (size_t i = 0; i < Rows; ++i)
    result[i][i] = diagonal[i];
  return result;
}

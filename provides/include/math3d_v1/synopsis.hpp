#pragma once

#include <cstdlib>

namespace math3d_v1 {

// trig.hpp ====================================================================

inline constexpr double pi = 3.141592653589793238462643383279502;

template <class Scalar> inline auto from_angle(Scalar angle);

// vec.hpp =====================================================================

template <class Scalar, size_t N> struct vec_t {
  using scalar_type = Scalar;

  Scalar values[N];

  auto dimensions() const;

  const auto &operator[](size_t i) const;
  auto &operator[](size_t i);
};

template <class BinOp, class SL, class SR, size_t N>
auto bin_op(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs);

template <class BinOp, class SL, class SR, size_t N>
auto bin_op(SL lhs, const vec_t<SR, N> &rhs);

template <class BinOp, class SL, class SR, size_t N>
auto bin_op(const vec_t<SL, N> &lhs, SR rhs);

#define MAKE(op)                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto op(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs);                   \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto op(SL lhs, const vec_t<SR, N> &rhs);                                    \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto op(const vec_t<SL, N> &lhs, SR rhs);

MAKE(operator+)
MAKE(operator-)
MAKE(operator*)
MAKE(operator/)
MAKE(operator%)
#undef MAKE

template <class SL, class SR, size_t N>
auto dot(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs);

template <class S, size_t N> auto norm(const vec_t<S, N> &v);

template <class S, size_t N> auto mag(const vec_t<S, N> &v);

template <class S, size_t N> auto zero_vec();

template <class SL, class SR, class ST, size_t N>
auto lerp(const vec_t<SL, N> &lhs, const vec_t<SR, N> &rhs, ST t);

template <class SL, class SR>
auto cross(const vec_t<SL, 3> &lhs, const vec_t<SR, 3> &rhs);

template <class S, size_t N> auto normalize(const vec_t<S, N> &v);

template <class S, size_t N> auto homogenize(const vec_t<S, N> &v);

template <size_t N, size_t I = 0, class S, size_t M>
auto sub(const vec_t<S, M> &v);

template <class Scalar, class... Scalars>
auto vec(Scalar value, Scalars... values);

template <class Scalar, size_t N, class... Scalars>
auto vec(const vec_t<Scalar, N> &v, Scalars... values);

// mtx.hpp =====================================================================

template <class Scalar, size_t R, size_t C = R> struct mtx_t {
  using scalar_type = Scalar;

  Scalar values[R][C];

  auto rows() const;
  auto columns() const;

  const auto &operator[](size_t i) const;
  auto &operator[](size_t i);
};

template <class BinOp, class SL, class SR, size_t R, size_t C>
auto bin_op(const mtx_t<SL, R, C> &lhs, const mtx_t<SR, R, C> &rhs);

#define MAKE(op)                                                               \
  template <class SL, class SR, size_t R, size_t C>                            \
  auto op(const mtx_t<SL, R, C> &lhs, const mtx_t<SR, R, C> &rhs);

MAKE(operator+)
MAKE(operator-)
#undef MAKE

template <class SL, class SR, size_t R, size_t RC, size_t C>
auto operator*(const mtx_t<SL, R, RC> &lhs, const mtx_t<SR, RC, C> &rhs);

template <class S, size_t R, size_t C> auto transpose(const mtx_t<S, R, C> &m);

template <class S, size_t N> auto set_identity(mtx_t<S, N> *out);

template <class S, size_t N> auto identity();

template <class S, size_t N, size_t M>
auto gauss_jordan(mtx_t<S, N> *a, mtx_t<S, N, M> *y);

template <class S, size_t N> auto inverse(const mtx_t<S, N> &m);

template <class SL, size_t N, size_t M, class SR>
auto operator*(const mtx_t<SL, N, M> &m, const vec_t<SR, M> &v);

template <size_t Rows, class... Scalars>
auto from_columns(const vec_t<Scalars, Rows> &...columns);

template <class Scalar, size_t Rows>
auto from_diagonal(const vec_t<Scalar, Rows> &diagonal);

// transform.hpp ===============================================================

template <class Scalar> auto perspective(Scalar fov, Scalar near, Scalar far);

template <class Scalar> auto translation(const vec_t<Scalar, 3> &trans);

template <class Scalar> auto scaling(const vec_t<Scalar, 3> &scale);

template <class Scalar> auto rotation(Scalar radians, size_t axis = 0);

} // namespace math3d_v1

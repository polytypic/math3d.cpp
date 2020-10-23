#pragma once

#include <cstdlib>

namespace math3d_v1 {

// =============================================================================

inline constexpr double pi = 3.141592653589793238462643383279502;

// vec.hpp =====================================================================

template <class Scalar, size_t N> struct vec {
  using scalar_type = Scalar;

  Scalar values[N];

  const Scalar &operator[](size_t i) const;
  Scalar &operator[](size_t i);
};

template <class BinOp, class SL, class SR, size_t N>
auto bin_op(const vec<SL, N> &lhs, const vec<SR, N> &rhs);

template <class BinOp, class SL, class SR, size_t N>
auto bin_op(SL lhs, const vec<SR, N> &rhs);

template <class BinOp, class SL, class SR, size_t N>
auto bin_op(const vec<SL, N> &lhs, SR rhs);

#define MAKE(name, op)                                                         \
  template <class SL, class SR, size_t N>                                      \
  auto op(const vec<SL, N> &lhs, const vec<SR, N> &rhs);                       \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto op(SL lhs, const vec<SR, N> &rhs);                                      \
                                                                               \
  template <class SL, class SR, size_t N>                                      \
  auto op(const vec<SL, N> &lhs, SR rhs);

MAKE(add_op, operator+)
MAKE(sub_op, operator-)
MAKE(mul_op, operator*)
MAKE(div_op, operator/)
MAKE(rem_op, operator%)
#undef MAKE

template <class SL, class SR, size_t N>
auto dot(const vec<SL, N> &lhs, const vec<SR, N> &rhs);

template <class S, size_t N> auto norm(const vec<S, N> &v);

template <class S, size_t N> auto mag(const vec<S, N> &v);

template <class S, size_t N> vec<S, N> zero_vec();

template <class SL, class SR, class ST, size_t N>
auto lerp(const vec<SL, N> &lhs, const vec<SR, N> &rhs, ST t);

template <class SL, class SR>
auto cross(const vec<SL, 3> &lhs, const vec<SR, 3> &rhs);

template <class S, size_t N> auto normalize(const vec<S, N> &v);

template <class S, size_t N> auto homogenize(const vec<S, N> &v);

template <size_t N, size_t I = 0, class S, size_t M>
auto sub(const vec<S, M> &v);

// mtx.hpp =====================================================================

template <class Scalar, size_t R, size_t C = R> struct mtx {
  using scalar_type = Scalar;

  Scalar values[R][C];

  const Scalar (&operator[](size_t i) const)[C];
  Scalar (&operator[](size_t i))[C];
};

template <class BinOp, class SL, class SR, size_t R, size_t C>
auto bin_op(const mtx<SL, R, C> &lhs, const mtx<SR, R, C> &rhs);

#define BIN_OP(name, op)                                                       \
  template <class SL, class SR, size_t R, size_t C>                            \
  auto op(const mtx<SL, R, C> &lhs, const mtx<SR, R, C> &rhs);

BIN_OP(add_op, operator+)
BIN_OP(sub_op, operator-)
#undef BIN_OP

template <class SL, class SR, size_t R, size_t RC, size_t C>
auto operator*(const mtx<SL, R, RC> &lhs, const mtx<SR, RC, C> &rhs);

template <class S, size_t R, size_t C> auto transpose(const mtx<S, R, C> &m);

template <class S, size_t N> auto set_identity(mtx<S, N> *out);

template <class S, size_t N> auto make_identity();

template <class S, size_t N, size_t M>
auto gauss_jordan(mtx<S, N> *a, mtx<S, N, M> *y);

template <class S, size_t N> auto inverse(const mtx<S, N> &m);

template <class SL, size_t N, size_t M, class SR>
auto operator*(const mtx<SL, N, M> &m, const vec<SR, M> &v);

// transform.hpp ===============================================================

template <class Scalar>
auto make_projection(Scalar fov, Scalar near, Scalar far);

template <class Scalar> auto make_translation(const vec<Scalar, 3> &trans);

template <class Scalar> auto make_scaling(const vec<Scalar, 3> &scale);

template <class Scalar> auto make_rotation(Scalar radians, size_t axis = 0);

} // namespace math3d_v1

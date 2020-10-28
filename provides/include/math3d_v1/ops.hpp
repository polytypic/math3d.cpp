#pragma once

namespace math3d_v1 {

#define MAKE(name, op)                                                         \
  struct name {                                                                \
    template <class L, class R>                                                \
    static inline auto apply(L lhs, R rhs) -> decltype(L() op R()) {           \
      return lhs op rhs;                                                       \
    }                                                                          \
  };

MAKE(add_op_t, +)
MAKE(sub_op_t, -)
MAKE(mul_op_t, *)
MAKE(div_op_t, /)
MAKE(rem_op_t, %)

MAKE(bxor_op_t, ^)
MAKE(band_op_t, &)
MAKE(bor_op_t, |)

MAKE(shl_op_t, <<)
MAKE(shr_op_t, >>)

MAKE(lt_op_t, <)
MAKE(lte_op_t, <=)
MAKE(eq_op_t, ==)
MAKE(neq_op_t, !=)
MAKE(gt_op_t, >)
MAKE(gte_op_t, >=)
#undef MAKE

#define MAKE(name, op)                                                         \
  struct name {                                                                \
    template <class T> static inline auto apply(T x) -> decltype(op T()) {     \
      return op x;                                                             \
    }                                                                          \
  };

MAKE(neg_op_t, -)
MAKE(bnot_op_t, ~)
MAKE(lnot_op_t, !)
#undef MAKE

} // namespace math3d_v1

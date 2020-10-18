#pragma once

namespace math3d_v1 {

#define MAKE(name, op)                                                         \
  struct name {                                                                \
    template <class L, class R>                                                \
    static inline auto apply(L lhs, R rhs) -> decltype(L() op R()) {           \
      return lhs op rhs;                                                       \
    }                                                                          \
  };

MAKE(add_op, +)
MAKE(sub_op, -)
MAKE(mul_op, *)
MAKE(div_op, /)
MAKE(rem_op, %)

MAKE(bxor_op, ^)
MAKE(band_op, &)
MAKE(bor_op, |)

MAKE(shl_op, <<)
MAKE(shr_op, >>)

MAKE(lt_op, <)
MAKE(lte_op, <=)
MAKE(eq_op, ==)
MAKE(neq_op, !=)
MAKE(gt_op, >)
MAKE(gte_op, >=)
#undef MAKE

#define MAKE(name, op)                                                         \
  struct name {                                                                \
    template <class T> static inline auto apply(T x) -> decltype(op T()) {     \
      return op x;                                                             \
    }                                                                          \
  };

MAKE(neg_op, -)
MAKE(bnot_op, ~)
MAKE(lnot_op, !)
#undef MAKE

} // namespace math3d_v1

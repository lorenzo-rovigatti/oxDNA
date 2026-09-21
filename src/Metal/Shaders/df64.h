/**
 * @file    df64.h
 * @brief   Double-float (df64) arithmetic for Metal shaders.
 *
 * A df64 number is an unevaluated sum hi+lo of two float32 values with
 * |lo| <= ulp(hi)/2, giving ~48 bits of effective precision (~2^-48
 * relative error). This is the Metal analogue of the hardware double
 * arithmetic used by the CUDA mixed-precision backend: Apple GPUs have no
 * `double` support in Metal Shading Language at all, so the `mixed`
 * precision tier keeps positions/velocities/angular momenta as df64 pairs
 * while forces/torques stay float32.
 *
 * Algorithms (error-free transformations):
 *  - TwoSum (Knuth 1969): exact rounded sum + rounding error.
 *  - QuickTwoSum: same when |a| >= |b| (used for renormalization).
 *  - TwoProd with FMA (Dekker 1971): exact rounded product + error.
 *    Metal provides `fma()` on float; Apple GPUs execute it at full speed.
 *  - df64+df64 addition (Dekker 1971, see also Thall 2006 "Extended-precision
 *    floating-point numbers for GPU computation").
 *
 * Only the operations needed by velocity-Verlet integration are provided:
 * df64+df64, df64+float, df64*float, float->df64 and df64->float.
 */

#ifndef DF64_H
#define DF64_H

#include <metal_stdlib>
using namespace metal;

/**
 * @brief A double-float number: value = hi + lo.
 */
struct df64 {
    float hi;
    float lo;

    df64() : hi(0.0f), lo(0.0f) {}
    df64(float h, float l) : hi(h), lo(l) {}
};

/// Construct a df64 from a single float (exact).
inline df64 df_from_f(float a) {
    df64 r;
    r.hi = a;
    r.lo = 0.0f;
    return r;
}

/// Round a df64 back to a single float (one rounding).
inline float df_to_f(df64 a) {
    return a.hi + a.lo;
}

/// Exact sum of two floats as (rounded sum, error). No ordering requirement.
inline df64 df_twosum(float a, float b) {
    float s = a + b;
    float v = s - a;
    float e = (a - (s - v)) + (b - v);
    df64 r;
    r.hi = s;
    r.lo = e;
    return r;
}

/// Renormalize a (sum, correction) pair with |sum| >= |correction|.
inline df64 df_renorm(float s, float e) {
    float hi = s + e;
    float lo = e - (hi - s);
    df64 r;
    r.hi = hi;
    r.lo = lo;
    return r;
}

/// Exact product of two floats as (rounded product, error) via FMA.
inline df64 df_twoprod(float a, float b) {
    float p = a * b;
    float e = fma(a, b, -p);
    df64 r;
    r.hi = p;
    r.lo = e;
    return r;
}

/// df64 + df64 (Dekker addition).
inline df64 df_add_df(df64 a, df64 b) {
    df64 s = df_twosum(a.hi, b.hi);
    float e = s.lo + a.lo + b.lo;
    return df_renorm(s.hi, e);
}

/// df64 + float.
inline df64 df_add_f(df64 a, float b) {
    df64 s = df_twosum(a.hi, b);
    float e = s.lo + a.lo;
    return df_renorm(s.hi, e);
}

/// df64 - float (via negation; exact).
inline df64 df_sub_f(df64 a, float b) {
    return df_add_f(a, -b);
}

/// df64 * float (TwoProd scaled + low-part correction, renormalized).
inline df64 df_mul_f(df64 a, float b) {
    df64 p = df_twoprod(a.hi, b);
    float e = p.lo + a.lo * b;
    return df_renorm(p.hi, e);
}

/// Split a double-precision host value into df64 (hi = float(v), lo = residual).
/// Used at init/sync time; the caller passes the float-cast value and the
/// residual computed on the host in double.
inline df64 df_split(float hi, float lo) {
    return df_renorm(hi, lo);
}

#endif /* DF64_H */

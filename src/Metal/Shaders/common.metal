/**
 * @file    common.metal
 * @brief   Common Metal shader functions and data structures
 *
 * Metal shader equivalent of common CUDA kernels
 */

#include <metal_stdlib>
using namespace metal;

#include "shader_utils.h"

// MetalBox, MetalBonds, MetalStressTensor defs removed (in shader_utils.h except StressTensor)
// MetalStressTensor was not in shader_utils.h? Let's check.
// I added MetalBox and MetalBonds.
// MetalStressTensor was in `metal_defs.h` (C++).
// common.metal had it too.
// I didn't verify MetalStressTensor in shader_utils.h.
// Let's check shader_utils.h content I wrote.
// I did NOT put MetalStressTensor in shader_utils.h.
// So I should keep MetalStressTensor here or move it.
// Ideally move it to shader_utils.h.

/**
 * @brief Stress tensor
 */
struct MetalStressTensor {
    m_number e[6];
};

/**
 * @brief Minimum image convention for periodic boundary conditions
 */
// minimum_image and distance_sqr moved to shader_utils.h

/**
 * @brief Safe arc cosine with clamping
 */
inline m_number safe_acos(m_number x) {
    if(x >= 1.0) return 0.0;
    if(x <= -1.0) return M_PI_F;
    return acos(x);
}

/**
 * @brief Cross product
 */
inline m_number3 cross_product(m_number3 a, m_number3 b) {
    return cross(a, b);
}

/**
 * @brief Normalize vector
 */
inline m_number3 normalize_vector(m_number3 v) {
    return normalize(v);
}

/**
 * @brief Atomic add for floating point (Metal doesn't have native atomic float add)
 */
inline void atomic_add_float(device atomic<float> *addr, float value) {
    // Use compare-and-swap loop
    float old_val = atomic_load_explicit(addr, memory_order_relaxed);
    float new_val;
    do {
        new_val = old_val + value;
    } while (!atomic_compare_exchange_weak_explicit(addr, &old_val, new_val,
                                                     memory_order_relaxed,
                                                     memory_order_relaxed));
}

/**
 * @brief Quaternion rotation
 */
inline m_number3 quaternion_rotate(m_number4 q, m_number3 v) {
    // q = (x, y, z, w) where w is the scalar part
    m_number3 qv = m_number3(q.x, q.y, q.z);
    m_number3 t = 2.0 * cross(qv, v);
    return v + q.w * t + cross(qv, t);
}

/**
 * @brief Parallel reduction for sum (threadgroup memory)
 */
template<typename T>
inline T threadgroup_reduce_add(T value, threadgroup T *shared_data, uint tid, uint threads) {
    shared_data[tid] = value;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Reduce within threadgroup
    for(uint s = threads / 2; s > 0; s >>= 1) {
        if(tid < s) {
            shared_data[tid] += shared_data[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    return shared_data[0];
}

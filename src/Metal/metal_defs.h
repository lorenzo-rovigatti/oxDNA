/*
 * metal_defs.h
 *
 *  Created for Metal GPU support on Apple Silicon
 *  Based on cuda_defs.h
 */

#ifndef SRC_METAL_METAL_DEFS_H_
#define SRC_METAL_METAL_DEFS_H_

#include <Metal/Metal.h>
#include <simd/simd.h>

/// Metal error checking macro
#define METAL_CHECK_ERROR(obj, msg)                                 \
  do {                                                              \
    if (!(obj)) {                                                   \
      fprintf(stderr, "Metal error at %s %d: %s\n",                \
             __FILE__, __LINE__, msg);                              \
      exit(EXIT_FAILURE);                                           \
    }                                                               \
  } while (0)

/// Thread execution width (equivalent to CUDA warp size, typically 32 for Apple GPUs)
#define METAL_SIMD_WIDTH 32

/// Useful macros for thread/threadgroup indexing
#define THREAD_INDEX_IN_THREADGROUP (thread_position_in_threadgroup)
#define THREADGROUP_INDEX (threadgroup_position_in_grid)
#define THREADS_PER_THREADGROUP (threads_per_threadgroup)
#define THREADGROUPS_IN_GRID (threadgroups_per_grid)

/// Global thread ID calculation
#define GLOBAL_THREAD_ID (threadgroup_position_in_grid * threads_per_threadgroup + thread_position_in_threadgroup)

/// Metal-specific arc cosine with clamping (similar to CUDA_LRACOS)
#define METAL_LRACOS(x) (((x) >= 1.0f) ? 0.0f : ((x) <= -1.0f) ? M_PI_F : acos(x))

/// Dot product macro for float3 vectors
#define METAL_DOT(a, b) (dot(a, b))

// Host-side scalar/vector types for the Metal backend.
//
// These are ALWAYS single precision: Apple GPUs have no hardware double
// precision (`double` is a compile error in Metal Shading Language), so the
// shader-side storage is float32 and the host mirrors must match it exactly.
// Higher accuracy is selected at run time with the `backend_precision` input
// key (float | mixed | hardmixed, see MD_MetalBackend.h), never at compile
// time. The legacy METAL_DOUBLE_PRECISION switch is therefore rejected.
#ifdef METAL_DOUBLE_PRECISION
#error "METAL_DOUBLE_PRECISION is not supported: Apple GPUs have no hardware double precision. Use backend_precision = 'mixed' (double-float emulation) or 'hardmixed' (CPU double integration) instead."
#endif
using m_number = float;
using m_number3 = simd_float3;
using m_number4 = simd_float4;
using m_quat = simd_float4;

/**
 * @brief Metal kernel configuration structure
 */
typedef struct Metal_kernel_cfg {
    uint32_t threadgroups_per_grid;
    uint32_t threads_per_threadgroup;
    uint32_t total_threads;
} Metal_kernel_cfg;

/**
 * @brief Neighbor bonds structure (3' and 5' directions)
 */
struct MetalBonds {
    int n3;
    int n5;
};

/**
 * @brief Edge bond structure for bond information
 */
struct MetalEdgeBond {
    int from;
    int to;
};

/**
 * @brief Stress tensor storage for Metal GPU
 */
struct MetalStressTensor {
    m_number e[6];

    MetalStressTensor() : e{0, 0, 0, 0, 0, 0} {}

    MetalStressTensor(m_number e0, m_number e1, m_number e2,
                      m_number e3, m_number e4, m_number e5) :
        e{e0, e1, e2, e3, e4, e5} {}

    inline MetalStressTensor operator+(const MetalStressTensor &other) const {
        return MetalStressTensor(
            e[0] + other.e[0],
            e[1] + other.e[1],
            e[2] + other.e[2],
            e[3] + other.e[3],
            e[4] + other.e[4],
            e[5] + other.e[5]);
    }

    inline void operator+=(const MetalStressTensor &other) {
        e[0] += other.e[0];
        e[1] += other.e[1];
        e[2] += other.e[2];
        e[3] += other.e[3];
        e[4] += other.e[4];
        e[5] += other.e[5];
    }

    inline void operator*=(const m_number other) {
        e[0] *= other;
        e[1] *= other;
        e[2] *= other;
        e[3] *= other;
        e[4] *= other;
        e[5] *= other;
    }

    inline void operator/=(const m_number other) {
        e[0] /= other;
        e[1] /= other;
        e[2] /= other;
        e[3] /= other;
        e[4] /= other;
        e[5] /= other;
    }
};

/**
 * @brief Metal buffer wrapper for type-safe buffer management
 */
template<typename T>
class MetalBuffer {
public:
    id<MTLBuffer> buffer;
    size_t size;
    size_t count;

    MetalBuffer() : buffer(nil), size(0), count(0) {}

    MetalBuffer(id<MTLDevice> device, size_t element_count) {
        count = element_count;
        size = element_count * sizeof(T);
        buffer = [device newBufferWithLength:size options:MTLResourceStorageModeShared];
        METAL_CHECK_ERROR(buffer, "Failed to allocate Metal buffer");
    }

    T* contents() {
        return static_cast<T*>([buffer contents]);
    }

    void release() {
        if (buffer) {
            buffer = nil;
        }
    }
};

#endif /* SRC_METAL_METAL_DEFS_H_ */

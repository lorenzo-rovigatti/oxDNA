/**
 * @file    shader_utils.h
 * @brief   Shared definitions for Metal shaders
 */

#ifndef SHADER_UTILS_H
#define SHADER_UTILS_H

#include <metal_stdlib>
using namespace metal;

// Precision selection.
//
// Apple GPUs have NO hardware double precision - `double` is a hard compile
// error in Metal Shading Language - so shader storage is always float32. Higher
// accuracy is provided at run time by the `backend_precision` input key:
//   * mixed     -> double-float (df64) arithmetic for positions/integration
//                  (see df64.metal)
//   * hardmixed -> the velocity-Verlet integration runs on the CPU in native
//                  double over the shared buffers
// The old compile-time -DMETAL_DOUBLE / METAL_DOUBLE_PRECISION switch only ever
// affected host-side types; it must never make the shader use `double`.
typedef float  m_number;
typedef float3 m_number3;
typedef float4 m_number4;

/**
 * @brief Simulation box structure
 */
struct MetalBox {
    m_number box_sides[3];
    m_number inv_sides[3];
};

/**
 * @brief Bond information
 */
struct MetalBonds {
    int n3;
    int n5;
};

/**
 * @brief One entry of the flat edge list (Metal_list = edge).
 */
struct MetalEdgeBond {
    int from;
    int to;
};

// Helper functions for Box
inline m_number3 minimum_image(m_number3 r, constant MetalBox &box) {
    m_number3 result = r;
    for(int i = 0; i < 3; i++) {
        result[i] -= box.box_sides[i] * rint(result[i] * box.inv_sides[i]);
    }
    return result;
}

inline m_number distance_sqr(m_number3 r1, m_number3 r2, constant MetalBox &box) {
    m_number3 dr = r1 - r2;
    dr = minimum_image(dr, box);
    return dot(dr, dr);
}

inline int compute_cell_index(constant int *N_cells_side, m_number4 r4, constant MetalBox &box) {
    m_number3 r = r4.xyz;
    int cx = (r.x * box.inv_sides[0] - floor(r.x * box.inv_sides[0])) * (1.0 - FLT_EPSILON) * N_cells_side[0];
    int cy = (r.y * box.inv_sides[1] - floor(r.y * box.inv_sides[1])) * (1.0 - FLT_EPSILON) * N_cells_side[1];
    int cz = (r.z * box.inv_sides[2] - floor(r.z * box.inv_sides[2])) * (1.0 - FLT_EPSILON) * N_cells_side[2];
    
    return (cz * N_cells_side[1] + cy) * N_cells_side[0] + cx;
}

inline int3 compute_cell_spl_idx(constant int *N_cells_side, m_number4 r4, constant MetalBox &box) {
    m_number3 r = r4.xyz;
    int cx = (r.x * box.inv_sides[0] - floor(r.x * box.inv_sides[0])) * (1.0 - FLT_EPSILON) * N_cells_side[0];
    int cy = (r.y * box.inv_sides[1] - floor(r.y * box.inv_sides[1])) * (1.0 - FLT_EPSILON) * N_cells_side[1];
    int cz = (r.z * box.inv_sides[2] - floor(r.z * box.inv_sides[2])) * (1.0 - FLT_EPSILON) * N_cells_side[2];
    
    return int3(cx, cy, cz);
}

#endif

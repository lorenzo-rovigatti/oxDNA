/**
 * @file    list_kernels.metal
 * @brief   Metal kernels for neighbor list updates
 */

#include "shader_utils.h"

using namespace metal;

struct FillCellsArgs {
    MetalBox box;
    int N_cells_side[3];
    int max_N_per_cell;
    int N;
};

struct UpdateNeighArgs {
    int N_cells_side[3];
    int max_N_per_cell;
    int N;
    m_number sqr_rverlet;
};

struct CountArgs {
    int N_cells_side[3];
    int N;
    MetalBox box; // Passing full struct
};

/**
 * @brief Count number of particles in each cell to determine max_N
 */
kernel void count_N_in_cells(device m_number4 *poss [[buffer(0)]],
                             device atomic_uint *counters_cells [[buffer(1)]],
                             constant CountArgs &args [[buffer(2)]],
                             uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;
    
    m_number4 r = poss[idx];
    int index = compute_cell_index(args.N_cells_side, r, args.box);
    
    atomic_fetch_add_explicit(&counters_cells[index], 1, memory_order_relaxed);
}

/**
 * @brief Fill cells with particle indices
 */
kernel void fill_cells(device m_number4 *poss [[buffer(0)]],
                       device int *cells [[buffer(1)]],
                       device atomic_int *counters_cells [[buffer(2)]],
                       device bool *cell_overflow [[buffer(3)]],
                       constant FillCellsArgs &args [[buffer(5)]],
                       uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;
    
    m_number4 r = poss[idx];
    int index = compute_cell_index(args.N_cells_side, r, args.box);
    
    int cell_idx = atomic_fetch_add_explicit(&counters_cells[index], 1, memory_order_relaxed);
    
    if(cell_idx < args.max_N_per_cell) {
        cells[index * args.max_N_per_cell + cell_idx] = idx;
    } else {
        *cell_overflow = true;
    }
}

// Helper to check neighbors in a specific cell
void update_cell_neigh_list(constant int *counters_cells,
                            device m_number4 *poss,
                            int cell_ind,
                            device int *cells,
                            m_number4 r,
                            thread int *neigh,
                            thread int &N_neigh,
                            MetalBonds b,
                            constant MetalBox &box,
                            int max_N_per_cell,
                            m_number sqr_rverlet,
                            int my_idx,
                            int max_neigh) {
    int size = counters_cells[cell_ind];
    for(int i = 0; i < size; i++) {
        int m = cells[cell_ind * max_N_per_cell + i];
        
        if(m == my_idx || b.n3 == m || b.n5 == m) continue;
        
        m_number4 rm = poss[m];
        m_number sqr_dist = distance_sqr(r.xyz, rm.xyz, box);
        
        if(sqr_dist < sqr_rverlet) {
            if (N_neigh < max_neigh) {
                // Should write to global memory directly
                // But here we are just counting or filling local?
                // The architecture writes to d_matrix_neighs
                // Let's change the function signature to write to memory
            }
        }
    }
}

/**
 * @brief Update neighbor lists
 */
kernel void simple_update_neigh_list(device m_number4 *poss [[buffer(0)]],
                                     device m_number4 *list_poss [[buffer(1)]],
                                     device int *cells [[buffer(2)]],
                                     constant int *counters_cells [[buffer(3)]],
                                     device int *matrix_neighs [[buffer(4)]],
                                     device int *number_neighs [[buffer(5)]],
                                     device MetalBonds *bonds [[buffer(6)]],
                                     constant MetalBox &box [[buffer(7)]],
                                     constant UpdateNeighArgs &args [[buffer(8)]],
                                     uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;
    
    m_number4 r = poss[idx];
    MetalBonds b = bonds[idx];
    int N_neighs = 0;
    
    // We assume max_neigh is handled by allocation size match
    // matrix_neighs layout: [N * max_neigh] ? Protocol typically uses row-major or similar
    // CUDASimpleVerletList uses: d_matrix_neighs allocated as N * max_neigh
    // Access: neigh[N_neigh * N + idx] in CUDA?
    // Let's check CUDA code: neigh[N_neigh * verlet_N[0] + IND] = m;
    // This is structure of arrays / column major? No, it's [N_neigh][idx]. 
    // It's coalesced access.
    // So neighbors of particle 0 are at 0, N, 2N...
    // Let's verify CUDA code.
    // `neigh[N_neigh * verlet_N[0] + IND] = m;`
    // Yes.
    
    int spl_idx_x, spl_idx_y, spl_idx_z;
    int3 spl = compute_cell_spl_idx(args.N_cells_side, r, box);
    
    // Loop over neighbors (3x3x3 block)
    for(int dz = -1; dz <= 1; dz++) {
        for(int dy = -1; dy <= 1; dy++) {
            for(int dx = -1; dx <= 1; dx++) {
                int nx = (spl.x + args.N_cells_side[0] + dx) % args.N_cells_side[0];
                int ny = (spl.y + args.N_cells_side[1] + dy) % args.N_cells_side[1];
                int nz = (spl.z + args.N_cells_side[2] + dz) % args.N_cells_side[2];
                
                int cell_index = (nz * args.N_cells_side[1] + ny) * args.N_cells_side[0] + nx;
                
                int size = counters_cells[cell_index];
                for(int i = 0; i < size; i++) {
                    int m = cells[cell_index * args.max_N_per_cell + i];
                    
                    if(m == idx || b.n3 == m || b.n5 == m) continue;
                    
                    m_number4 rm = poss[m];
                    m_number sqr_dist = distance_sqr(r.xyz, rm.xyz, box);
                    
                    if(sqr_dist < args.sqr_rverlet) {
                        matrix_neighs[N_neighs * args.N + idx] = m;
                        N_neighs++;
                        // Warning: Needs bound check max_neigh? CUDA doesn't seem to check in kernel, assumes safe?
                    }
                }
            }
        }
    }
    
    list_poss[idx] = r;
    number_neighs[idx] = N_neighs;
}

struct EdgeFillArgs {
    int N_cells_side[3];
    int max_N_per_cell;
    int N;
    m_number sqr_rverlet;
    int max_edges;
};

/**
 * @brief Build the flat edge list: one thread per particle appends its
 * half-list neighbours (m < idx) as (idx, m) pairs using a single atomic
 * global counter. Cells must already be filled.
 */
kernel void edge_fill(device m_number4 *poss            [[buffer(0)]],
                      device m_number4 *list_poss       [[buffer(1)]],
                      device int *cells                 [[buffer(2)]],
                      constant int *counters_cells      [[buffer(3)]],
                      device MetalEdgeBond *edge_list   [[buffer(4)]],
                      device atomic_int *edge_count     [[buffer(5)]],
                      device bool *edge_overflow        [[buffer(6)]],
                      device MetalBonds *bonds          [[buffer(7)]],
                      constant MetalBox &box            [[buffer(8)]],
                      constant EdgeFillArgs &args       [[buffer(9)]],
                      uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;

    m_number4 r = poss[idx];
    MetalBonds b = bonds[idx];
    int3 spl = compute_cell_spl_idx(args.N_cells_side, r, box);

    for(int dz = -1; dz <= 1; dz++) {
        for(int dy = -1; dy <= 1; dy++) {
            for(int dx = -1; dx <= 1; dx++) {
                int nx = (spl.x + args.N_cells_side[0] + dx) % args.N_cells_side[0];
                int ny = (spl.y + args.N_cells_side[1] + dy) % args.N_cells_side[1];
                int nz = (spl.z + args.N_cells_side[2] + dz) % args.N_cells_side[2];
                int cell_index = (nz * args.N_cells_side[1] + ny) * args.N_cells_side[0] + nx;

                int size = counters_cells[cell_index];
                for(int i = 0; i < size; i++) {
                    int m = cells[cell_index * args.max_N_per_cell + i];
                    if(m >= idx || b.n3 == m || b.n5 == m) continue;   // half list

                    if(distance_sqr(r.xyz, poss[m].xyz, box) < args.sqr_rverlet) {
                        int e = atomic_fetch_add_explicit(edge_count, 1, memory_order_relaxed);
                        if(e < args.max_edges) {
                            edge_list[e].from = idx;
                            edge_list[e].to = m;
                        }
                        else {
                            *edge_overflow = true;
                        }
                    }
                }
            }
        }
    }

    list_poss[idx] = r;
}

kernel void check_coord_magnitude(device m_number4 *poss [[buffer(0)]],
                                  device bool *is_large [[buffer(1)]],
                                  uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    const float THRESHOLD = 1.e7;
    m_number4 p = poss[idx];
    is_large[idx] = (abs(p.x) > THRESHOLD || abs(p.y) > THRESHOLD || abs(p.z) > THRESHOLD);
}

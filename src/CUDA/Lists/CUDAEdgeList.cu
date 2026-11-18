/*
 * CUDAEdgeList.cu
 *
 * Author: Subhajit Roy
 *
 * GPU direct edge-list implementation.
 *
 * Builds unique (IND, m) pairs with IND > m via a two-pass approach:
 *   Pass 1: count_edges   — record per-particle edge count
 *   scan  : exclusive_scan — compute per-particle offsets
 *   Pass 2: fill_edges    — write edge_bond entries at computed offsets
 *
 * Memory savings vs CUDASimpleVerletList+use_edge:
 *   No d_matrix_neighs (N * max_neigh * sizeof(int)) allocation.
 *   d_edge_list size = N_edges * sizeof(edge_bond) ~= N * avg_neigh/2 * 8 bytes.
 */

#include "CUDAEdgeList.h"
#include "CUDA_edge.cuh"
#include "../../Utilities/oxDNAException.h"
#include "../../Utilities/Utils.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Particles/BaseParticle.h"
#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../CUDAUtils.h"

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/extrema.h>
#include <thrust/transform.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>

CUDAEdgeList::CUDAEdgeList() {
    _use_edge = true;  // always produces edge list
    _cells_kernel_cfg.threads_per_block = 0;
    _N_cells = _old_N_cells = N_edges = -1;
}

CUDAEdgeList::~CUDAEdgeList() {
}

void CUDAEdgeList::clean() {
    if(_counters_cells_tex != 0) {
        cudaDestroyTextureObject(_counters_cells_tex);
        _counters_cells_tex = 0;
    }
    if(_d_cells != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_cells));
        CUDA_SAFE_CALL(cudaFree(_d_counters_cells));
        CUDA_SAFE_CALL(cudaFreeHost(_d_cell_overflow));
        _d_cells = _d_counters_cells = nullptr;
    }
    if(_d_edge_counts != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_edge_counts));
        CUDA_SAFE_CALL(cudaFree(_d_edge_offsets));
        _d_edge_counts = _d_edge_offsets = nullptr;
    }
    if(d_edge_list != nullptr) {
        CUDA_SAFE_CALL(cudaFree(d_edge_list));
        d_edge_list = nullptr;
    }
}

void CUDAEdgeList::get_settings(input_file &inp) {
    getInputBool(&inp, "cells_auto_optimisation", &_auto_optimisation, 0);
    getInputNumber(&inp, "verlet_skin", &_verlet_skin, 1);
    getInputNumber(&inp, "max_density_multiplier", &_max_density_multiplier, 0);

    // CUDAEdgeList only ever populates d_edge_list/N_edges: it never
    // allocates the per-particle d_matrix_neighs/d_number_neighs that the
    // non-edge force kernels need. Without use_edge = true those kernels
    // silently fall back to their O(N^2) all-pairs path (still correct,
    // just needlessly slow), so require use_edge explicitly instead.
    bool use_edge = false;
    getInputBool(&inp, "use_edge", &use_edge, 0);
    if(!use_edge) {
        throw oxDNAException("CUDA_list = edge requires use_edge = true");
    }
}

// --------------------------------------------------------------------------
// cell layout helpers (mirror CUDASimpleVerletList internals)
// --------------------------------------------------------------------------

__global__ void _edge_count_N_in_cells(c_number4 *poss, uint *counters_cells,
                                         int N_cells_side[3], int N, CUDABox box) {
    if(IND >= N) return;
    c_number4 r = poss[IND];
    int index = box.compute_cell_index(N_cells_side, r);
    atomicInc((uint *) &counters_cells[index], N);
}

void CUDAEdgeList::_compute_N_cells_side(int N_cells_side[3], c_number min_cell_size) {
    c_number4 box_sides_n4 = _h_cuda_box->box_sides();
    c_number box_sides[3] = { box_sides_n4.x, box_sides_n4.y, box_sides_n4.z };
    c_number max_factor = pow(2. * _N / _h_cuda_box->V(), 1. / 3.);

    for(int i = 0; i < 3; i++) {
        N_cells_side[i] = (int) (floor(box_sides[i] / min_cell_size) + 0.1);
        if(N_cells_side[i] < 3) N_cells_side[i] = 3;
        if(_auto_optimisation && N_cells_side[i] > ceil(max_factor * box_sides[i])) {
            N_cells_side[i] = ceil(max_factor * box_sides[i]);
        }
    }
}

int CUDAEdgeList::_largest_N_in_cells(c_number4 *poss, c_number min_cell_size) {
    int N = CONFIG_INFO->N();

    int *N_cells_side;
    CUDA_SAFE_CALL(cudaMallocHost(&N_cells_side, sizeof(int) * 3, cudaHostAllocDefault));
    _compute_N_cells_side(N_cells_side, min_cell_size);
    int N_cells = N_cells_side[0] * N_cells_side[1] * N_cells_side[2];

    uint *counters_cells;
    CUDA_SAFE_CALL(cudaMalloc(&counters_cells, (size_t) N_cells * sizeof(uint)));
    CUDA_SAFE_CALL(cudaMemset(counters_cells, 0, N_cells * sizeof(uint)));

    int tpb = 64;
    int blocks = N / tpb + ((N % tpb == 0) ? 0 : 1);
    _edge_count_N_in_cells<<<blocks, tpb>>>(poss, counters_cells, N_cells_side, N, *_h_cuda_box);
    CUT_CHECK_ERROR("_edge_count_N_in_cells error");

    thrust::device_ptr<uint> dev_ptr(counters_cells);
    int max_N = *thrust::max_element(dev_ptr, dev_ptr + N_cells);

    CUDA_SAFE_CALL(cudaFreeHost(N_cells_side));
    CUDA_SAFE_CALL(cudaFree(counters_cells));

    return max_N;
}

void CUDAEdgeList::_init_cells(c_number4 *poss) {
    _compute_N_cells_side(_N_cells_side, std::sqrt(_sqr_rverlet));
    _N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

    if(_old_N_cells != -1 && _N_cells != _old_N_cells) {
        CUDA_SAFE_CALL(cudaFree(_d_cells));
        CUDA_SAFE_CALL(cudaFree(_d_counters_cells));
        _d_cells = _d_counters_cells = nullptr;
        cudaDestroyTextureObject(_counters_cells_tex);
        OX_DEBUG("CUDAEdgeList: re-allocating cells, from %d to %d\n", _old_N_cells, _N_cells);
    }

    if(_d_cells == nullptr) {
        bool deallocate = false;
        if(poss == nullptr) {
            deallocate = true;
            int N = CONFIG_INFO->N();
            std::vector<c_number4> host_positions;
            host_positions.reserve(N);
            for(auto p : CONFIG_INFO->particles()) {
                c_number4 pos({(c_number) p->pos[0], (c_number) p->pos[1], (c_number) p->pos[2], 0.});
                host_positions.push_back(pos);
            }
            CUDA_SAFE_CALL(cudaMalloc(&poss, (size_t) N * sizeof(c_number4)));
            CUDA_SAFE_CALL(cudaMemcpy(poss, host_positions.data(), sizeof(c_number4) * N, cudaMemcpyHostToDevice));
        }

        _max_N_per_cell = std::round(_max_density_multiplier * _largest_N_in_cells(poss, std::sqrt(_sqr_rverlet)));
        if(_max_N_per_cell > _N) _max_N_per_cell = _N + 1;
        if(_max_N_per_cell < 5) _max_N_per_cell = 5;

        CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_counters_cells, (size_t) _N_cells * sizeof(int)));
        CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_cells, (size_t) _N_cells * _max_N_per_cell * sizeof(int)));

        CUDA_SAFE_CALL(cudaMemcpyToSymbol(edge_N_cells_side,    _N_cells_side,     3 * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(edge_max_N_per_cell, &_max_N_per_cell,  sizeof(int)));

        if(_counters_cells_tex != 0) cudaDestroyTextureObject(_counters_cells_tex);
        GpuUtils::init_texture_object(&_counters_cells_tex,
            cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindSigned),
            _d_counters_cells, _N_cells);

        if(deallocate) CUDA_SAFE_CALL(cudaFree(poss));
    }

    _old_N_cells = _N_cells;
}

void CUDAEdgeList::_realloc_edge_list(int new_capacity) {
    if(d_edge_list != nullptr) CUDA_SAFE_CALL(cudaFree(d_edge_list));
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&d_edge_list, (size_t) new_capacity * sizeof(edge_bond)));
    _edge_list_capacity = new_capacity;
}

void CUDAEdgeList::init(int N, c_number rcut, CUDABox *h_cuda_box, CUDABox *d_cuda_box) {
    CUDABaseList::init(N, rcut, h_cuda_box, d_cuda_box);

    c_number rverlet   = rcut + 2 * _verlet_skin;
    _sqr_rverlet       = SQR(rverlet);
    _sqr_verlet_skin   = SQR(_verlet_skin);

    _init_cells();

    // estimate initial edge list capacity: N * max_neigh / 2
    int est_max_neigh = std::min((int)(4 * M_PI * _max_N_per_cell / 3.), N - 1);
    int initial_capacity = std::max(N * est_max_neigh / 2, N);
    _realloc_edge_list(initial_capacity);

    // per-particle edge counts and offsets (N+1 for the scan)
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_edge_counts,  (size_t)(_N + 1) * sizeof(int)));
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_edge_offsets, (size_t)(_N + 1) * sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(edge_N, &_N, sizeof(int)));

    CUDA_SAFE_CALL(cudaMallocHost(&_d_cell_overflow, sizeof(bool), cudaHostAllocDefault));
    _d_cell_overflow[0] = false;

    if(_cells_kernel_cfg.threads_per_block == 0) _cells_kernel_cfg.threads_per_block = 64;
    _cells_kernel_cfg.blocks.x = _N / _cells_kernel_cfg.threads_per_block +
                                  ((_N % _cells_kernel_cfg.threads_per_block == 0) ? 0 : 1);
    _cells_kernel_cfg.blocks.y = _cells_kernel_cfg.blocks.z = 1;

    float f_copy = _sqr_rverlet;
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(edge_sqr_rverlet, &f_copy, sizeof(float)));

    OX_LOG(Logger::LOG_INFO, "CUDAEdgeList: rverlet=%.3g, initial edge capacity=%d, max_N_per_cell=%d",
           (double)rverlet, initial_capacity, _max_N_per_cell);
}

void CUDAEdgeList::update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds) {
    _init_cells(poss);
    CUDA_SAFE_CALL(cudaMemset(_d_counters_cells, 0, _N_cells * sizeof(int)));

    // fill cells
    edge_fill_cells
        <<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
        (poss, _d_cells, _d_counters_cells, _d_cell_overflow, _d_cuda_box);
    CUT_CHECK_ERROR("edge_fill_cells error");
    cudaDeviceSynchronize();

    if(_d_cell_overflow[0]) {
        throw oxDNAException("CUDAEdgeList: cell overflow (max_N_per_cell=%d). "
                              "Try max_density_multiplier = 10 in the input file.", _max_N_per_cell);
    }

    // pass 1: count edges per particle
    count_edges
        <<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
        (_counters_cells_tex, poss, list_poss, _d_cells, _d_edge_counts, bonds, _d_cuda_box);
    CUT_CHECK_ERROR("count_edges error");

    // exclusive scan: d_edge_offsets[i] = sum of d_edge_counts[0..i-1]
    // write the extra element so offsets[N] = total edges
    thrust::device_ptr<int> d_counts(_d_edge_counts);
    thrust::device_ptr<int> d_offsets(_d_edge_offsets);
    // store N_edges at offsets[N]: pad d_edge_counts[N]=0 then scan N+1 elements
    CUDA_SAFE_CALL(cudaMemset(_d_edge_counts + _N, 0, sizeof(int)));
    thrust::exclusive_scan(d_counts, d_counts + _N + 1, d_offsets);

    // read total edge count from device
    int total_edges = 0;
    CUDA_SAFE_CALL(cudaMemcpy(&total_edges, _d_edge_offsets + _N, sizeof(int), cudaMemcpyDeviceToHost));
    N_edges = total_edges;

    // grow edge list buffer if needed (add 10% headroom)
    if(total_edges > _edge_list_capacity) {
        int new_cap = (int)(total_edges * 1.1) + 64;
        OX_DEBUG("CUDAEdgeList: growing edge list from %d to %d", _edge_list_capacity, new_cap);
        _realloc_edge_list(new_cap);
    }

    // pass 2: fill edge list
    fill_edges
        <<<_cells_kernel_cfg.blocks, _cells_kernel_cfg.threads_per_block>>>
        (_counters_cells_tex, poss, _d_cells, _d_edge_offsets, d_edge_list, bonds, _d_cuda_box);
    CUT_CHECK_ERROR("fill_edges error");
}

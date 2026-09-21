/**
 * @file    MetalEdgeList.h
 * @brief   Flat edge (pair) list for the Metal backend.
 *
 * Reuses the cell machinery of MetalSimpleVerletList but, instead of a
 * per-particle neighbour matrix, builds one contiguous array of unique
 * (from, to) pairs (from > to). The DNA force kernel then runs one thread per
 * edge and evaluates each non-bonded pair once, atomically distributing the
 * force/torque to both partners - the Metal analogue of CUDA's use_edge /
 * CUDAEdgeList path. Selected with `Metal_list = edge`.
 */

#ifndef METALEDGELIST_H_
#define METALEDGELIST_H_

#include "MetalSimpleVerletList.h"

class MetalEdgeList : public MetalSimpleVerletList {
protected:
    id<MTLComputePipelineState> _edge_fill_pso = nil;

    id<MTLBuffer> _d_edge_overflow = nil;
    int _max_edges = 0;

    void _alloc_edge_buffers(int max_edges);

public:
    MetalEdgeList() = default;
    virtual ~MetalEdgeList();

    void metal_init(int N, m_number rcut, MetalBox *h_metal_box, id<MTLBuffer> d_metal_box,
                    id<MTLDevice> device, id<MTLLibrary> library) override;
    void update(id<MTLBuffer> poss, id<MTLBuffer> list_poss, id<MTLBuffer> bonds) override;
    void clean() override;

    bool is_edge_list() override { return true; }
};

#endif /* METALEDGELIST_H_ */

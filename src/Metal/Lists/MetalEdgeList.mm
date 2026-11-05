/**
 * @file MetalEdgeList.mm
 */

#include "MetalEdgeList.h"
#include "../../Utilities/oxDNAException.h"

MetalEdgeList::~MetalEdgeList() {
    clean();
}

void MetalEdgeList::clean() {
    MetalSimpleVerletList::clean();
    _edge_fill_pso = nil;
    _d_edge_overflow = nil;
    d_edge_list = nil;
    d_n_edges = nil;
}

void MetalEdgeList::_alloc_edge_buffers(int max_edges) {
    _max_edges = max_edges;
    d_edge_list = MetalUtils::allocate_buffer<MetalEdgeBond>(_device, _max_edges, MTLResourceStorageModePrivate);
    OX_LOG(Logger::LOG_INFO, "Metal edge list: capacity %d edges (%.1f MB)",
           _max_edges, _max_edges * (double) sizeof(MetalEdgeBond) / 1048576.);
}

void MetalEdgeList::metal_init(int N, m_number rcut, MetalBox *h_metal_box, id<MTLBuffer> d_metal_box,
                               id<MTLDevice> device, id<MTLLibrary> library) {
    MetalSimpleVerletList::metal_init(N, rcut, h_metal_box, d_metal_box, device, library);

    // The edge path does not use the per-particle neighbour matrix.
    d_matrix_neighs = nil;
    d_number_neighs = nil;

    NSError *error = nil;
    _edge_fill_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"edge_fill"] error:&error];
    if(!_edge_fill_pso) {
        throw oxDNAException("Failed to create edge_fill pipeline: %s", [[error localizedDescription] UTF8String]);
    }

    d_n_edges = MetalUtils::allocate_buffer<int>(_device, 1, MTLResourceStorageModeShared);
    _d_edge_overflow = MetalUtils::allocate_buffer<bool>(_device, 1, MTLResourceStorageModeShared);

    // Initial capacity estimate; grows on overflow.
    int est = std::min(_max_neigh, 96);
    _alloc_edge_buffers((int) ((long long) N * std::max(est, 8)));
}

void MetalEdgeList::update(id<MTLBuffer> poss, id<MTLBuffer> list_poss, id<MTLBuffer> bonds) {
    _run_fill_cells(poss);

    struct EdgeFillArgs {
        int N_cells_side[3];
        int max_N_per_cell;
        int N;
        m_number sqr_rverlet;
        int max_edges;
    } args;
    args.N_cells_side[0] = _N_cells_side[0];
    args.N_cells_side[1] = _N_cells_side[1];
    args.N_cells_side[2] = _N_cells_side[2];
    args.max_N_per_cell = _max_N_per_cell;
    args.N = _N;
    args.sqr_rverlet = _sqr_rverlet;

    for(int attempt = 0; attempt < 4; attempt++) {
        args.max_edges = _max_edges;

        int zero = 0;
        MetalUtils::copy_to_device(d_n_edges, &zero, 1);
        bool of = false;
        MetalUtils::copy_to_device(_d_edge_overflow, &of, 1);

        id<MTLCommandQueue> q = [_device newCommandQueue];
        id<MTLCommandBuffer> cb = [q commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];

        [enc setComputePipelineState:_edge_fill_pso];
        [enc setBuffer:poss offset:0 atIndex:0];
        [enc setBuffer:list_poss offset:0 atIndex:1];
        [enc setBuffer:_d_cells offset:0 atIndex:2];
        [enc setBuffer:_d_counters_cells offset:0 atIndex:3];
        [enc setBuffer:d_edge_list offset:0 atIndex:4];
        [enc setBuffer:d_n_edges offset:0 atIndex:5];
        [enc setBuffer:_d_edge_overflow offset:0 atIndex:6];
        [enc setBuffer:bonds offset:0 atIndex:7];
        [enc setBuffer:_d_metal_box offset:0 atIndex:8];
        [enc setBytes:&args length:sizeof(args) atIndex:9];

        [enc dispatchThreadgroups:MTLSizeMake(_cells_kernel_cfg.threadgroups_per_grid, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(_cells_kernel_cfg.threads_per_threadgroup, 1, 1)];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];

        MetalUtils::copy_from_device(&of, _d_edge_overflow, 1);
        MetalUtils::copy_from_device(&N_edges, d_n_edges, 1);

        if(!of) {
            return;
        }
        // grow and retry
        int needed = (int) (N_edges * 1.3) + 1024;
        _alloc_edge_buffers(std::max(needed, _max_edges * 2));
    }
    throw oxDNAException("Metal edge list kept overflowing (last count %d)", N_edges);
}

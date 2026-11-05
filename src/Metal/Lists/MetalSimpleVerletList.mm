/*
 * MetalSimpleVerletList.mm
 *
 *  Created for Metal backend
 */

#include "MetalSimpleVerletList.h"
#include "../../Utilities/oxDNAException.h"
#include "../../Utilities/Utils.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Particles/BaseParticle.h"
// #include "../MetalUtils.h" // Already included via MetalBaseList

#include <cmath>
#include <vector>
#include <algorithm>

MetalSimpleVerletList::MetalSimpleVerletList() {
	_cells_kernel_cfg.threads_per_threadgroup = 0;
	_use_edge = false;
	_N_cells = _old_N_cells = N_edges = -1;
    _device = nil;
    _library = nil;
}

MetalSimpleVerletList::~MetalSimpleVerletList() {
    clean();
}

void MetalSimpleVerletList::clean() {
	if(_d_cells != nil) {
        _d_cells = nil;
        _d_counters_cells = nil;
        d_matrix_neighs = nil;
        d_number_neighs = nil;
        _d_cell_overflow = nil;
	}

	if(_use_edge && d_edge_list != nil) {
        d_edge_list = nil;
        _d_number_neighs_no_doubles = nil;
	}
}

void MetalSimpleVerletList::get_settings(input_file &inp) {
	getInputBool(&inp, "cells_auto_optimisation", &_auto_optimisation, 0);
	getInputBool(&inp, "print_problematic_ids", &_print_problematic_ids, 0);
	getInputNumber(&inp, "verlet_skin", &_verlet_skin, 1);
	getInputFloat(&inp, "max_density_multiplier", &_max_density_multiplier, 0);
    OX_LOG(Logger::LOG_INFO, "MetalList: max_density_multiplier set to %f", _max_density_multiplier);
	bool want_edge = false;
	getInputBool(&inp, "use_edge", &want_edge, 0);
	if(want_edge) {
		// use_edge is a CUDA neighbour-list layout optimisation. The Metal DNA
		// kernel uses the plain neighbour matrix and the edge kernels are not
		// implemented, so ignore the key instead of crashing on a nil pipeline.
		OX_LOG(Logger::LOG_WARNING, "use_edge is not supported by the Metal backend and will be ignored");
	}
	_use_edge = false;
}

int MetalSimpleVerletList::_largest_N_in_cells(id<MTLBuffer> poss, m_number min_cell_size) {
    if(!poss) return 0;
    
	int N = CONFIG_INFO->N();

	int N_cells_side[3];
	_compute_N_cells_side(N_cells_side, min_cell_size);
	int N_cells = N_cells_side[0] * N_cells_side[1] * N_cells_side[2];

    // Allocate temp buffer for counters
    id<MTLBuffer> counters_cells = MetalUtils::allocate_buffer<uint>( _device, N_cells, MTLResourceStorageModeShared);
    memset(counters_cells.contents, 0, N_cells * sizeof(uint));

	int threads_per_threadgroup = 64;
	int threadgroups = N / threads_per_threadgroup + ((N % threads_per_threadgroup == 0) ? 0 : 1);

    id<MTLCommandQueue> commandQueue = [_device newCommandQueue];
    id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
    
    [computeEncoder setComputePipelineState:_count_N_in_cells_pso];
    [computeEncoder setBuffer:poss offset:0 atIndex:0];
    [computeEncoder setBuffer:counters_cells offset:0 atIndex:1];
    
    // Pass N_cells_side, N, and box data
    MetalBox::BoxData box_data = _h_metal_box->get_box_data();
    struct {
        int N_cells_side[3];
        int N;
        MetalBox::BoxData box;
    } args;
    args.N_cells_side[0] = N_cells_side[0];
    args.N_cells_side[1] = N_cells_side[1];
    args.N_cells_side[2] = N_cells_side[2];
    args.N = N;
    args.box = box_data;
    
    [computeEncoder setBytes:&args length:sizeof(args) atIndex:2];

    [computeEncoder dispatchThreadgroups:MTLSizeMake(threadgroups, 1, 1) threadsPerThreadgroup:MTLSizeMake(threads_per_threadgroup, 1, 1)];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    uint *host_counters = (uint*)[counters_cells contents];
    uint max_N = 0;
    for(int i=0; i<N_cells; i++) {
        if(host_counters[i] > max_N) max_N = host_counters[i];
    }
    
    // Release temp buffer implicitly handled by ARC if local variable? No, I allocated with MetalUtils? 
    // MetalUtils uses [device newBuffer...]. ARC handles it.
    
	return (int)max_N;
}

void MetalSimpleVerletList::_compute_N_cells_side(int N_cells_side[3], m_number min_cell_size) {
    if(!_h_metal_box) return;
    
	LR_vector box_sides = {_h_metal_box->get_box_data().box_sides[0], _h_metal_box->get_box_data().box_sides[1], _h_metal_box->get_box_data().box_sides[2]};
	// c_number max_factor = pow(2. * _N / _h_metal_box->V(), 1. / 3.);
    // Accessing V() via host helper
    m_number vol = _h_metal_box->V();
    m_number max_factor = pow(2. * _N / vol, 1. / 3.);

	for(int i = 0; i < 3; i++) {
        m_number side = (i==0 ? box_sides.x : (i==1 ? box_sides.y : box_sides.z));
        
		N_cells_side[i] = (int) (floor(side / min_cell_size) + 0.1);
		if(N_cells_side[i] < 3) {
			N_cells_side[i] = 3;
		}
		if(_auto_optimisation && N_cells_side[i] > ceil(max_factor * side)) {
			N_cells_side[i] = ceil(max_factor * side);
		}
	}
}

void MetalSimpleVerletList::_init_cells(id<MTLBuffer> poss) {
	_compute_N_cells_side(_N_cells_side, std::sqrt(_sqr_rverlet));
	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	if(_old_N_cells != -1 && _N_cells != _old_N_cells) {
        _d_cells = nil;
        _d_counters_cells = nil;
		OX_DEBUG("Re-allocating cells on GPU, from %d to %d\n", _old_N_cells, _N_cells);
	}

	if(_d_cells == nil) {
		bool deallocate = false;
		if(poss == nil) {
			deallocate = true;
			int N = CONFIG_INFO->N();
            // Allocate temporary buffer if poss is nil
            // But we need data to compute largest N.
            // If poss is nil, we might need to copy from host or something.
            // CUDA version copies from host CONFIG_INFO if poss is null.
            // We can do the same.
            
            size_t size = N * sizeof(m_number4);
            poss = MetalUtils::allocate_buffer<m_number4>(_device, N);
            m_number4 *host_ptr = (m_number4*) poss.contents;
			for(auto p : CONFIG_INFO->particles()) {
                host_ptr[p->index] = (m_number4){(m_number) p->pos[0], (m_number) p->pos[1], (m_number) p->pos[2], 0.};
			}
		}

		_max_N_per_cell = std::round(_max_density_multiplier * _largest_N_in_cells(poss, std::sqrt(_sqr_rverlet)));
		if(_max_N_per_cell > _N) {
			_max_N_per_cell = _N;
		}
		if(_max_N_per_cell < 5) {
			_max_N_per_cell = 5;
		}

        _d_counters_cells = MetalUtils::allocate_buffer<int>(_device, _N_cells, MTLResourceStorageModeShared);
        _d_cells = MetalUtils::allocate_buffer<int>(_device, _N_cells * _max_N_per_cell, MTLResourceStorageModePrivate);
        
		if(deallocate) {
            poss = nil; // ARC handles release
		}
	}

	_old_N_cells = _N_cells;
}

void MetalSimpleVerletList::metal_init(int N, m_number rcut, MetalBox *h_metal_box, id<MTLBuffer> d_metal_box, id<MTLDevice> device, id<MTLLibrary> library) {
	MetalBaseList::metal_init(N, rcut, h_metal_box, d_metal_box, device, library);
    _device = device;
    _library = library;
    
    // Create Pipeline States
    NSError *error = nil;
    
    _fill_cells_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"fill_cells"] error:&error];
    if (!_fill_cells_pso) printf("Error creating fill_cells PSO: %s\n", [[error localizedDescription] UTF8String]);
    
    _update_neigh_list_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"simple_update_neigh_list"] error:&error];
    if (!_update_neigh_list_pso) printf("Error creating simple_update_neigh_list PSO: %s\n", [[error localizedDescription] UTF8String]);

    _count_N_in_cells_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"count_N_in_cells"] error:&error];
    if (!_count_N_in_cells_pso) printf("Error creating count_N_in_cells PSO: %s\n", [[error localizedDescription] UTF8String]);
    
    if(_use_edge) {
         _edge_update_neigh_list_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"edge_update_neigh_list"] error:&error];
         _compress_matrix_neighs_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"compress_matrix_neighs"] error:&error];
    }
    
    _check_coord_magnitude_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"check_coord_magnitude"] error:&error];

	m_number rverlet = rcut + 2 * _verlet_skin;
	_sqr_rverlet = SQR(rverlet);
	_sqr_verlet_skin = SQR(_verlet_skin);
	_vec_size = N * sizeof(m_number4);

	_init_cells();

	OX_LOG(Logger::LOG_INFO, "Metal Cells mem: %.2lf MBs, lists mem: %.2lf MBs", (double) _N_cells*(1 + _max_N_per_cell) * sizeof(int)/1048576., (double) _N * (1 + _max_neigh) * sizeof(int)/1048576.);

	_max_neigh = std::min((int) (4 * M_PI * _max_N_per_cell / 3.), N - 1);
	OX_LOG(Logger::LOG_INFO, "Metal max_neigh: %d, max_N_per_cell: %d, N_cells: %d (per side: %d %d %d)", _max_neigh, _max_N_per_cell, _N_cells, _N_cells_side[0], _N_cells_side[1], _N_cells_side[2]);

    d_number_neighs = MetalUtils::allocate_buffer<int>(_device, _N, MTLResourceStorageModeShared); // Should be private ideally but maybe we read it back?
    d_matrix_neighs = MetalUtils::allocate_buffer<int>(_device, _N * _max_neigh, MTLResourceStorageModePrivate);

    _d_cell_overflow = MetalUtils::allocate_buffer<bool>(_device, 1, MTLResourceStorageModeShared);
    bool fail = false;
    MetalUtils::copy_to_device(_d_cell_overflow, &fail, 1);

	if(_use_edge) {
        d_edge_list = MetalUtils::allocate_buffer<MetalEdgeBond>(_device, _N * _max_neigh, MTLResourceStorageModePrivate);
        _d_number_neighs_no_doubles = MetalUtils::allocate_buffer<int>(_device, _N + 1, MTLResourceStorageModeShared);
	}

	if(_cells_kernel_cfg.threads_per_threadgroup == 0) _cells_kernel_cfg.threads_per_threadgroup = 64;
	_cells_kernel_cfg.threadgroups_per_grid = _N / _cells_kernel_cfg.threads_per_threadgroup + ((_N % _cells_kernel_cfg.threads_per_threadgroup == 0) ? 0 : 1);
    _cells_kernel_cfg.total_threads = _N;

	OX_DEBUG("Cells kernel cfg: threads_per_threadgroup = %d, threadgroups = %d", _cells_kernel_cfg.threads_per_threadgroup, _cells_kernel_cfg.threadgroups_per_grid);
}

std::vector<int> MetalSimpleVerletList::is_large(id<MTLBuffer> data) {
    // Implementing check for problematic particles
    // For now, simpler implementation or skip
    return std::vector<int>();
}

bool MetalSimpleVerletList::lists_are_old(id<MTLBuffer> poss, id<MTLBuffer> list_poss) {
    // Shared-storage buffers on Apple Silicon: read directly on the host.
    if(poss == nil || list_poss == nil) return true;

    const m_number4 *p = (const m_number4 *) poss.contents;
    const m_number4 *lp = (const m_number4 *) list_poss.contents;
    if(p == nullptr || lp == nullptr) return true;

    const m_number sqr_skin = _sqr_verlet_skin;
    for(int i = 0; i < _N; i++) {
        m_number dx = p[i].x - lp[i].x;
        m_number dy = p[i].y - lp[i].y;
        m_number dz = p[i].z - lp[i].z;
        if(dx * dx + dy * dy + dz * dz > sqr_skin) {
            return true;
        }
    }
    return false;
}

// Assigns every particle to its cell (fill_cells kernel) and throws on a cell
// overflow. Shared by the verlet and edge list update paths.
void MetalSimpleVerletList::_run_fill_cells(id<MTLBuffer> poss) {
    _init_cells(poss); // Check if cells need resize

    bool fail = false;
    MetalUtils::copy_to_device(_d_cell_overflow, &fail, 1);
    memset(_d_counters_cells.contents, 0, _N_cells * sizeof(int));

    id<MTLCommandQueue> commandQueue = [_device newCommandQueue];
    id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];

    [computeEncoder setComputePipelineState:_fill_cells_pso];
    [computeEncoder setBuffer:poss offset:0 atIndex:0];
    [computeEncoder setBuffer:_d_cells offset:0 atIndex:1];
    [computeEncoder setBuffer:_d_counters_cells offset:0 atIndex:2];
    [computeEncoder setBuffer:_d_cell_overflow offset:0 atIndex:3];

    struct {
        MetalBox::BoxData box;
        int N_cells_side[3];
        int max_N_per_cell;
        int N;
    } args;
    args.box = _h_metal_box->get_box_data();
    args.N_cells_side[0] = _N_cells_side[0];
    args.N_cells_side[1] = _N_cells_side[1];
    args.N_cells_side[2] = _N_cells_side[2];
    args.max_N_per_cell = _max_N_per_cell;
    args.N = _N;
    [computeEncoder setBytes:&args length:sizeof(args) atIndex:5];

    [computeEncoder dispatchThreadgroups:MTLSizeMake(_cells_kernel_cfg.threadgroups_per_grid, 1, 1) threadsPerThreadgroup:MTLSizeMake(_cells_kernel_cfg.threads_per_threadgroup, 1, 1)];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    MetalUtils::copy_from_device(&fail, _d_cell_overflow, 1);
    if(fail) {
        throw oxDNAException("A cell contains more than _max_n_per_cell (%d) particles:", _max_N_per_cell);
    }
}

void MetalSimpleVerletList::update(id<MTLBuffer> poss, id<MTLBuffer> list_poss, id<MTLBuffer> bonds) {
    _run_fill_cells(poss);

    id<MTLCommandQueue> commandQueue = [_device newCommandQueue];
    id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
    
	if(_use_edge) {
		// edge_update_neigh_list
        [computeEncoder setComputePipelineState:_edge_update_neigh_list_pso];
        // Set args similar to simple... plus edge specific
        // ... implementation pending
	}
	else {
		// simple_update_neigh_list
        [computeEncoder setComputePipelineState:_update_neigh_list_pso];
        [computeEncoder setBuffer:poss offset:0 atIndex:0];
        [computeEncoder setBuffer:list_poss offset:0 atIndex:1];
        [computeEncoder setBuffer:_d_cells offset:0 atIndex:2];
        [computeEncoder setBuffer:_d_counters_cells offset:0 atIndex:3]; // Using buffer instead of texture
        [computeEncoder setBuffer:d_matrix_neighs offset:0 atIndex:4];
        [computeEncoder setBuffer:d_number_neighs offset:0 atIndex:5];
        [computeEncoder setBuffer:bonds offset:0 atIndex:6];
        [computeEncoder setBuffer:_d_metal_box offset:0 atIndex:7];
        
        struct {
            int N_cells_side[3];
            int max_N_per_cell;
            int N;
            m_number sqr_rverlet;
        } args;
        args.N_cells_side[0] = _N_cells_side[0];
        args.N_cells_side[1] = _N_cells_side[1];
        args.N_cells_side[2] = _N_cells_side[2];
        args.max_N_per_cell = _max_N_per_cell;
        args.N = _N;
        args.sqr_rverlet = _sqr_rverlet;
        
        [computeEncoder setBytes:&args length:sizeof(args) atIndex:8];
        
        [computeEncoder dispatchThreadgroups:MTLSizeMake(_cells_kernel_cfg.threadgroups_per_grid, 1, 1) threadsPerThreadgroup:MTLSizeMake(_cells_kernel_cfg.threads_per_threadgroup, 1, 1)];
	}
    
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

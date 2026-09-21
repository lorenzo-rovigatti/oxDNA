/*
 * MetalDNAInteraction.mm
 *
 *  Created for Metal backend
 */

#include "MetalDNAInteraction.h"
#include "MetalCPUForceFallback.h"
#include "../Lists/MetalSimpleVerletList.h"
#include "../Lists/MetalNoList.h" 

// Define parameter struct to match shader
struct DNAInteractionParams {
    float F1_EPS[50];
    float F1_SHIFT[50];
    float F1_A[2];
    float F1_RC[2];
    float F1_R0[2];
    float F1_BLOW[2];
    float F1_BHIGH[2];
    float F1_RLOW[2];
    float F1_RHIGH[2];
    float F1_RCLOW[2];
    float F1_RCHIGH[2];
    
    float F2_K[2];
    float F2_RC[2];
    float F2_R0[2];
    float F2_BLOW[2];
    float F2_BHIGH[2];
    float F2_RLOW[2];
    float F2_RHIGH[2];
    float F2_RCLOW[2];
    float F2_RCHIGH[2];
    
    float F5_PHI_A[4];
    float F5_PHI_B[4];
    float F5_PHI_XC[4];
    float F5_PHI_XS[4];
    
    float hb_multiplier;
    float T;
    
    // Debye Huckel
    float dh_RC;
    float dh_RHIGH;
    float dh_prefactor;
    float dh_B;
    float dh_minus_kappa;
    int dh_half_charged_ends;
    
    // Flags
    int grooving;
    int use_oxDNA2_coaxial_stacking;
    int use_oxDNA2_FENE;
    int use_debye_huckel;
    int use_mbf;
    float mbf_xmax;
    float mbf_finf;
    float sqr_rcut;
};

MetalDNAInteraction::MetalDNAInteraction() {
	_edge_compatible = true;
    _d_is_strand_end = nil;
    _d_dna_params = nil;
    _d_init_args = nil;
    _command_queue = nil;
}

MetalDNAInteraction::~MetalDNAInteraction() {
    _d_is_strand_end = nil;
    _d_dna_params = nil;
}

void MetalDNAInteraction::get_settings(input_file &inp) {
    // Metal specific settings if any
    DNAInteraction::get_settings(inp);
    
    // Handle oxDNA2 logic similar to CUDA/CPU
	std::string inter_type;
	if(getInputString(&inp, "interaction_type", inter_type, 0) == KEY_FOUND) {
		if(inter_type.compare("DNA2") == 0) {
            // Logic copied from CUDADNAInteraction/DNA2Interaction
			_use_debye_huckel = true;
			_use_oxDNA2_coaxial_stacking = true;
			_use_oxDNA2_FENE = true;
            F2_K[1] = CXST_K_OXDNA2;
			_debye_huckel_half_charged_ends = true;
			this->_grooving = true;
            
			getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
			getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);

			_debye_huckel_lambdafactor = 0.3616455f;
			getInputFloat(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0);

			_debye_huckel_prefactor = 0.0543f;
			getInputFloat(&inp, "dh_strength", &_debye_huckel_prefactor, 0);
		}
	}
}

void MetalDNAInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
    DNAInteraction::check_input_sanity(particles);
}

void MetalDNAInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
    DNAInteraction::read_topology(N_strands, particles);
}

void MetalDNAInteraction::metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library) {
    MetalBaseInteraction::metal_init(N, device, library);
    DNAInteraction::init();
    
    // Initialize parameters buffer
    _d_dna_params = MetalUtils::allocate_buffer<DNAInteractionParams>(_device, 1, MTLResourceStorageModeShared);
    DNAInteractionParams *params = (DNAInteractionParams*)[_d_dna_params contents];
    
    params->hb_multiplier = this->_hb_multiplier;
    params->T = this->_T;
    
#define COPY_ARR(DEST, SRC, N) for(int i=0; i<N; i++) DEST[i] = SRC[i];
#define COPY_MAT(DEST, SRC, D1, D2, D3) \
    for(int i=0; i<D1; i++) for(int j=0; j<D2; j++) for(int k=0; k<D3; k++) \
        DEST[i*D2*D3 + j*D3 + k] = SRC[i][j][k];

    COPY_MAT(params->F1_EPS, this->F1_EPS, 2, 5, 5);
    COPY_MAT(params->F1_SHIFT, this->F1_SHIFT, 2, 5, 5);
    
    COPY_ARR(params->F1_A, this->F1_A, 2);
    COPY_ARR(params->F1_RC, this->F1_RC, 2);
    COPY_ARR(params->F1_R0, this->F1_R0, 2);
    COPY_ARR(params->F1_BLOW, this->F1_BLOW, 2);
    COPY_ARR(params->F1_BHIGH, this->F1_BHIGH, 2);
    COPY_ARR(params->F1_RLOW, this->F1_RLOW, 2);
    COPY_ARR(params->F1_RHIGH, this->F1_RHIGH, 2);
    COPY_ARR(params->F1_RCLOW, this->F1_RCLOW, 2);
    COPY_ARR(params->F1_RCHIGH, this->F1_RCHIGH, 2);
    
    COPY_ARR(params->F2_K, this->F2_K, 2);
    COPY_ARR(params->F2_RC, this->F2_RC, 2);
    COPY_ARR(params->F2_R0, this->F2_R0, 2);
    COPY_ARR(params->F2_BLOW, this->F2_BLOW, 2);
    COPY_ARR(params->F2_BHIGH, this->F2_BHIGH, 2);
    COPY_ARR(params->F2_RLOW, this->F2_RLOW, 2);
    COPY_ARR(params->F2_RHIGH, this->F2_RHIGH, 2);
    COPY_ARR(params->F2_RCLOW, this->F2_RCLOW, 2);
    COPY_ARR(params->F2_RCHIGH, this->F2_RCHIGH, 2);
    
    COPY_ARR(params->F5_PHI_A, this->F5_PHI_A, 4);
    COPY_ARR(params->F5_PHI_B, this->F5_PHI_B, 4);
    COPY_ARR(params->F5_PHI_XC, this->F5_PHI_XC, 4);
    COPY_ARR(params->F5_PHI_XS, this->F5_PHI_XS, 4);
    
    // Debye-Huckel initialization logic (simplified port from CUDA)
    if (_use_debye_huckel) {
		m_number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
		_debye_huckel_RHIGH = 3.0 * lambda;
		_minus_kappa = -1.0 / lambda;

		m_number x = _debye_huckel_RHIGH;
		m_number q = _debye_huckel_prefactor;
		m_number l = lambda;
		_debye_huckel_B = -(exp(-x / l) * q * q * (x + l) * (x + l)) / (-4. * x * x * x * l * l * q);
		_debye_huckel_RC = x * (q * x + 3. * q * l) / (q * (x + l));

		m_number debyecut;
		if(this->_grooving) {
			debyecut = 2.0f * sqrt(SQR(POS_MM_BACK1) + SQR(POS_MM_BACK2)) + _debye_huckel_RC;
		}
		else {
			debyecut = 2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
		}
		if(debyecut > this->_rcut) {
			this->_rcut = debyecut;
			this->_sqr_rcut = debyecut * debyecut;
		}
    }
    
    params->dh_RC = _debye_huckel_RC;
    params->dh_RHIGH = _debye_huckel_RHIGH;
    params->dh_prefactor = _debye_huckel_prefactor;
    params->dh_B = _debye_huckel_B;
    params->dh_minus_kappa = _minus_kappa;
    params->dh_half_charged_ends = _debye_huckel_half_charged_ends ? 1 : 0;
    
    params->grooving = _grooving ? 1 : 0;
    params->use_oxDNA2_coaxial_stacking = _use_oxDNA2_coaxial_stacking ? 1 : 0;
    params->use_oxDNA2_FENE = this->_use_oxDNA2_FENE ? 1 : 0;
    params->use_debye_huckel = _use_debye_huckel ? 1 : 0;

    // Max backbone force ("relax"-style FENE cap). The kernel needs xmax/finf,
    // matching CUDA's _bonded_part.
    params->use_mbf = this->_use_mbf ? 1 : 0;
    params->mbf_xmax = this->_mbf_xmax;
    params->mbf_finf = this->_mbf_finf;

    // COM cutoff^2 (after the Debye-Huckel block may have grown _rcut). The
    // kernel uses this to skip Verlet-list entries that sit in the skin shell.
    params->sqr_rcut = (float) this->_sqr_rcut;
    
    // PSOs
    NSError *error = nil;
    _dna_forces_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"dna_forces"] error:&error];
    if(!_dna_forces_pso) {
        throw oxDNAException("Failed creating dna_forces pipeline: %s", [[error localizedDescription] UTF8String]);
    }

    _dna_forces_edge_nonbonded_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"dna_forces_edge_nonbonded"] error:&error];
    if(!_dna_forces_edge_nonbonded_pso) {
        throw oxDNAException("Failed creating dna_forces_edge_nonbonded pipeline: %s", [[error localizedDescription] UTF8String]);
    }
    _dna_forces_edge_bonded_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"dna_forces_edge_bonded"] error:&error];
    if(!_dna_forces_edge_bonded_pso) {
        throw oxDNAException("Failed creating dna_forces_edge_bonded pipeline: %s", [[error localizedDescription] UTF8String]);
    }

    _init_DNA_strand_ends_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"init_DNA_strand_ends"] error:&error];
    if(!_init_DNA_strand_ends_pso) {
        throw oxDNAException("Failed creating init_DNA_strand_ends pipeline: %s", [[error localizedDescription] UTF8String]);
    }
    
    // Allocate init args (N)
    struct InitStrandArgs { int N; };
    InitStrandArgs args;
    args.N = _N;
    _d_init_args = MetalUtils::allocate_buffer<InitStrandArgs>(_device, 1);
    MetalUtils::copy_to_device(_d_init_args, &args, 1);
}

void MetalDNAInteraction::_on_T_update() {
    metal_init(_N, _device, _library);
}

void MetalDNAInteraction::_init_strand_ends(id<MTLBuffer> d_bonds) {
    if(!_d_is_strand_end) {
        _d_is_strand_end = MetalUtils::allocate_buffer<int>(_device, _N, MTLResourceStorageModePrivate);
    }
    
    id<MTLCommandQueue> queue = [_device newCommandQueue];
    id<MTLCommandBuffer> commandBuffer = [queue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
    
    [computeEncoder setComputePipelineState:_init_DNA_strand_ends_pso];
    [computeEncoder setBuffer:_d_is_strand_end offset:0 atIndex:0];
    [computeEncoder setBuffer:d_bonds offset:0 atIndex:1];
    struct { int N; } args = { _N };
    [computeEncoder setBytes:&args length:sizeof(args) atIndex:2];
    
    int tpb = 64;
    int blocks = _N / tpb + ((_N % tpb == 0) ? 0 : 1);
    [computeEncoder dispatchThreadgroups:MTLSizeMake(blocks, 1, 1) threadsPerThreadgroup:MTLSizeMake(tpb, 1, 1)];
    
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalDNAInteraction::compute_forces(MetalBaseList *lists, id<MTLBuffer> d_poss, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_bonds, id<MTLBuffer> d_box, id<MTLBuffer> d_energies) {
    if(_use_cpu_fallback) {
        MetalCPUForceFallback::compute(_N, d_poss, d_orientations, d_forces, d_torques, d_energies);
        return;
    }

    if(!_d_is_strand_end) _init_strand_ends(d_bonds);

    // Check if pipelines are ready
    if(!_dna_forces_pso) {
        throw oxDNAException("MetalDNAInteraction non-bonded pipeline not initialized");
    }

    @autoreleasepool {
        if(!_command_queue) {
            _command_queue = [_device newCommandQueue];
        }
        id<MTLCommandBuffer> commandBuffer = [_command_queue commandBuffer];
        this->process_dna_force_kernel(commandBuffer, lists, d_poss, d_orientations, d_forces, d_torques, d_bonds, d_box, d_energies);
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
    }
}

void MetalDNAInteraction::encode_forces(id<MTLCommandBuffer> command_buffer, MetalBaseList *lists, id<MTLBuffer> d_poss,
                                        id<MTLBuffer> d_orientations, id<MTLBuffer> forces, id<MTLBuffer> torques,
                                        id<MTLBuffer> bonds, id<MTLBuffer> metal_box, id<MTLBuffer> energies) {
    if(_use_cpu_fallback) {
        MetalCPUForceFallback::compute(_N, d_poss, d_orientations, forces, torques, energies);
        return;
    }
    if(!_d_is_strand_end) _init_strand_ends(bonds);
    if(!_dna_forces_pso) {
        throw oxDNAException("MetalDNAInteraction non-bonded pipeline not initialized");
    }
    this->process_dna_force_kernel(command_buffer, lists, d_poss, d_orientations, forces, torques, bonds, metal_box, energies);
}

void MetalDNAInteraction::process_dna_force_kernel(id<MTLCommandBuffer> commandBuffer,
                                                   MetalBaseList *list,
                                                   id<MTLBuffer> poss,
                                                   id<MTLBuffer> orientations,
                                                   id<MTLBuffer> forces,
                                                   id<MTLBuffer> torques,
                                                   id<MTLBuffer> bonds,
                                                   id<MTLBuffer> metal_box,
                                                   id<MTLBuffer> energies) {
    const MTLSize tpb = MTLSizeMake(128, 1, 1); // 128 measured fastest (register-bound)

    @autoreleasepool {
        if(list->is_edge_list()) {
            // one thread per edge: non-bonded pair computed once, atomically
            // scattered to both partners (lab-frame torque)
            id<MTLComputeCommandEncoder> nb = [commandBuffer computeCommandEncoder];
            [nb setComputePipelineState:_dna_forces_edge_nonbonded_pso];
            [nb setBuffer:poss offset:0 atIndex:0];
            [nb setBuffer:orientations offset:0 atIndex:1];
            [nb setBuffer:forces offset:0 atIndex:2];
            [nb setBuffer:torques offset:0 atIndex:3];
            [nb setBuffer:list->d_edge_list offset:0 atIndex:4];
            [nb setBuffer:list->d_n_edges offset:0 atIndex:5];
            [nb setBuffer:bonds offset:0 atIndex:6];
            [nb setBuffer:_d_dna_params offset:0 atIndex:7];
            [nb setBuffer:metal_box offset:0 atIndex:8];
            int n_edges = (list->N_edges > 0) ? list->N_edges : 1;
            [nb dispatchThreads:MTLSizeMake(n_edges, 1, 1) threadsPerThreadgroup:tpb];
            [nb endEncoding];

            // one thread per particle: add bonded terms, transform torque to body frame
            id<MTLComputeCommandEncoder> bd = [commandBuffer computeCommandEncoder];
            [bd setComputePipelineState:_dna_forces_edge_bonded_pso];
            [bd setBuffer:poss offset:0 atIndex:0];
            [bd setBuffer:orientations offset:0 atIndex:1];
            [bd setBuffer:forces offset:0 atIndex:2];
            [bd setBuffer:torques offset:0 atIndex:3];
            [bd setBuffer:bonds offset:0 atIndex:4];
            [bd setBuffer:_d_dna_params offset:0 atIndex:5];
            [bd setBuffer:metal_box offset:0 atIndex:6];
            [bd setBuffer:_d_init_args offset:0 atIndex:7];
            [bd dispatchThreads:MTLSizeMake(_N, 1, 1) threadsPerThreadgroup:tpb];
            [bd endEncoding];
            return;
        }

        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:_dna_forces_pso];

        [encoder setBuffer:poss offset:0 atIndex:0];
        [encoder setBuffer:orientations offset:0 atIndex:1];
        [encoder setBuffer:forces offset:0 atIndex:2];
        [encoder setBuffer:torques offset:0 atIndex:3];
        [encoder setBuffer:list->d_matrix_neighs offset:0 atIndex:4];
        [encoder setBuffer:list->d_number_neighs offset:0 atIndex:5];
        [encoder setBuffer:bonds offset:0 atIndex:6];
        [encoder setBuffer:_d_dna_params offset:0 atIndex:7];
        [encoder setBuffer:metal_box offset:0 atIndex:8];
        [encoder setBuffer:_d_init_args offset:0 atIndex:9];

        // Energy buffer at index 10
        if(energies) {
             [encoder setBuffer:energies offset:0 atIndex:10];
        }

        [encoder dispatchThreads:MTLSizeMake(_N, 1, 1) threadsPerThreadgroup:tpb];
        [encoder endEncoding];
    }
}

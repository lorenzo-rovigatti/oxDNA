/*
 * CUDAPHBInteraction.h
 *
 *  Created on: ongoing
 *      Author: subho
 */

#ifndef CUDAPHBINTERACTION_H_
#define CUDAPHBINTERACTION_H_

#define CUDA_MAX_PATCHES 30 // subho: Be careful, increase this if needed

#include "CUDABaseInteraction.h"
#include "../../Interactions/PHBInteraction.h"
#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

/**
 * @Subhajit-Roy-Partho
 * @brief This class implements both patchy and the helix bundles interactions in CUDA
 */

__constant__ float rcut2;  // cut-off distance squared
__constant__ int exclusionType; // 0 for Linear, 1 for Cubic, 2 for Hard
__constant__ float sigma; // 
__constant__ float patchyB; // Controls the stiffness of exe volume and in case of hard the power over (sigma/r).
__constant__ int NumPatches[3]; // Total number of patches for particle type 0 is ico or main particle, 1 is helix and 2 is no patches
__constant__ float4 basePatchConfig[3][CUDA_MAX_PATCHES]; // Same as patchConfig
__constant__ float patchyRcutSqr;
__constant__ float patchyAlpha;
__constant__ float patchyEpsilon;
__constant__ float hardVolCutoff;
__constant__ float connections[MAX_Particle][MAX_NEIGHBOURS];

class CUDAPHBInteraction: public CUDABaseInteraction, public PHBInteraction {
public:
    CUDAPHBInteraction();
    virtual ~CUDAPHBInteraction();
    void get_settings(input_file &inp) override;
    void cuda_init(int N) override;
    c_number get_cuda_rcut() {
		return this->get_rcut();
	}
    void compute_forces(CUDABaseList *lists,c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torque, LR_bonds *d_bonds, CUDABox *d_box);

};

// void CUDAinterParticleInteraction();

#endif // CUDAPHBINTERACTION_H_
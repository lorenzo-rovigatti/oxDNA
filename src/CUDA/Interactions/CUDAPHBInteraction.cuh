/*
 * CUDAPHBInteraction.h
 *
 *  Created on: ongoing
 *      Author: subho
 */

#ifndef CUDAPHBINTERACTION_H_
#define CUDAPHBINTERACTION_H_

#define CUDA_MAX_PATCHES 50 // subho: Be careful, increase this if needed
#define MAX_Particle 1000
#define MAX_Neighbour 50

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
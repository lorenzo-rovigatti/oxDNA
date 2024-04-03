#include "CUDACCGInteration.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

__device__ float GPUk[MAXparticles][MAXneighbour+1];//spring constant
__device__ float GPUro[MAXparticles][MAXneighbour+1];//rarius of the spring
__device__ float GPUconnection[MAXparticles][MAXneighbour+1];//first intger will state number of connections
__device__ float GPUtopology[MAXparticles][4];

__constant__ float GPUpatchyEpsilon;
__constant__ float GPUpatchySigma;
__constant__ float GPUpatchyRstar;
__constant__ float GPUpatchyB;
__constant__ float GPUpatchyRc;
__constant__ float GPUpatchyRcut;
__constant__ float GPUpatchyAlpha;
__constant__ float GPUpatchyEcutoff;

__constant__ int MD_N[1];


__device__ void repulsive(c_number prefactor,c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc){
    c_number rnorm = CUDA_DOT(r,r);
    F.x = F.y = F.z = F.w = (c_number) 0.f;
    if(rnorm < SQR(rc)){
        if(rnorm > SQR(rstar)){
            c_number rmod = sqrt(rnorm);
            c_number rrc = rmod -rc;
            c_number fmod = 2.f*prefactor*b*rrc/rmod;
            F.x = r.x * fmod;
            F.y = r.y * fmod;
            F.z = r.z * fmod;
            F.w = prefactor*b*SQR(rrc);
            return;
        }
        c_number lj_part = CUB(SQR(sigma)/rnorm);
        c_number fmod = 24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / rnorm;
        F.x = r.x * fmod;
		F.y = r.y * fmod;
		F.z = r.z * fmod;
        F.w = 4.f * prefactor * (SQR(lj_part) - lj_part);
        return;
    }
}

__device__ void CCGforce(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *number_neighs, CUDABox *box){
    if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

    int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
    for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND) {
			c_number4 qpos = poss[k_index];

			GPU_quat qo = orientations[k_index];
			get_vectors_from_quat(qo, b1, b2, b3);
			// _patchy_two_body_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, bonds, k_index, box);
		}
	}

    forces[IND] = F;
    torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);

}


CUDACCGInteraction::CUDACCGInteraction(){

};

CUDACCGInteraction::~CUDACCGInteraction(){

};

void CUDACCGInteraction::get_settings(input_file &inp){
    CCGInteraction::get_settings(inp);
};

void CUDACCGInteraction::cuda_init(int N){ // N is the number of particles
    CUDABaseInteraction::cuda_init(N); //Firing up cuda 
    CCGInteraction::init();
    std::vector<BaseParticle *> particles(N); // Initiate the base particles
	int Nstrands;
    CCGInteraction::read_topology(&Nstrands,particles);
    size_t k_size = N * sizeof(c_number);
};

void CUDACCGInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_qorientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box){

};
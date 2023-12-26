#include "CUDACCGInteration.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

CUDACCGInteraction::CUDACCGInteraction(){

};

CUDACCGInteraction::~CUDACCGInteraction(){

};

void CUDACCGInteraction::get_settings(input_file &inp){

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
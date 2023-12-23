#include "CUDACCGInteration.h"

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
    if(rnorm < SQR(rc)){
        if(rnorm > SQR(rstar)){
            rmod = 
        }
    }
}
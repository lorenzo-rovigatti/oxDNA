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
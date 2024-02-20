#include "CUDAPHBInteraction.cuh"


CUDAPHBInteraction::CUDAPHBInteraction() {
    _edge_compatible = true;
}

CUDAPHBInteraction::~CUDAPHBInteraction() {
}

void CUDAPHBInteraction::get_settings(input_file &inp) {
    PHBInteraction::get_settings(inp);
}


void CUDAPHBInteraction::cuda_init(int N){
    CUDABaseInteraction::cuda_init(N);
    PHBInteraction::init();
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(rcut2, &_sqr_rcut, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(sigma, &patchySigma, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyB, &patchyB, sizeof(float)));
    
}

void CUDAPHBInteraction::compute_forces(CUDABaseList *lists,c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torque, LR_bonds *d_bonds, CUDABox *d_box) {
    // std::cout<<"rcut"<<_rcut<<std::endl;
    return;
}

__device__ void CUDAinterParticleInteraction(c_number4 &pPos, c_number4 &qPos, c_number4 &a1, c_number4 &a2, c_number4 &a3, c_number4 &b1, c_number4 &b2, c_number4 &b3, c_number4 &f, c_number4 &tor, CUDABox *box ){ //__device__ can only be called from a __global__ function, this only run on the GPU
    c_number4 r = box->minimum_image(pPos,qPos); //minimum distance in periodic boundary conditions
    c_number r2 = CUDA_DOT(r, r); // r^2
    if(r2 >= rcut2) return; //if r^2 is greater than total cutoff stop



}

__device__ void CUDAexeVolLin(c_number prefactor,c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc){
    auto r2 = CUDA_DOT(r,r);
    if(r2 <SQR(rc)){
        if(r2 > SQR(rstar)){
            c_number rmod = sqrt(r2);
            c_number rrc = rmod -rc;
            c_number fmod = 2.f*prefactor*b*rrc/rmod;
            F.x -= r.x * fmod;
            F.y -= r.y * fmod;
            F.z -= r.z * fmod;
            F.w += prefactor*b*SQR(rrc);
            return;
        }
        c_number lj_part = CUB(SQR(sigma)/r2);
        c_number fmod = 24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / r2;
        F.x -= r.x * fmod;
        F.y -= r.y * fmod;
        F.z -= r.z * fmod;
        F.w += 4.f * prefactor * (SQR(lj_part) - lj_part);
        return;
    }
}

__device__ void CUDAexeVolCub(c_number prefactor,c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc){
    auto r2 = CUDA_DOT(r,r);
    if(r2 <SQR(rc)){
        if(r2 > SQR(rstar)){
            c_number rmod = sqrt(r2);
            c_number rrc = rmod -rc;
            c_number fmod = 2.f*prefactor*b*CUB(rrc)/rmod;
            F.x -= r.x * fmod;
            F.y -= r.y * fmod;
            F.z -= r.z * fmod;
            F.w += prefactor*b*SQR(SQR(rrc));
            return;
        }
        c_number lj_part = CUB(SQR(sigma)/r2);
        c_number fmod = 24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / r2;
        F.x -= r.x * fmod;
        F.y -= r.y * fmod;
        F.z -= r.z * fmod;
        F.w += 4.f * prefactor * (SQR(lj_part) - lj_part);
        return;
    }
}

__device__ void CUDAexeVolHard(c_number4 &r,c_number4 &F){
    c_number r2= CUDA_DOT(r,r);
    c_number part = powf(sigma/r2, patchyB*0.5);
    F.w += part;
    c_number fmod= patchyB*part/r2;
    F.x -= r.x*fmod;
    F.y -= r.y*fmod;
    F.z -= r.z*fmod;
}

__device__ void CUDAspring(){
    return;
}

__device__ void CUDAbondedTwist(){
    return;
}

__device__ void CUDAbondedDoubleBending(){
    return;
}

__device__ void CUDAbondedAlignment(){
    return;
}

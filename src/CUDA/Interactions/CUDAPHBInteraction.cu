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
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(NumPatches, &NumPatches, sizeof(int)*3));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(basePatchConfig, &basePatchConfig, sizeof(float4)*3*CUDA_MAX_PATCHES));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyRcutSqr, &patchyRcutSqr, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyAlpha, &patchyAlpha, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyEpsilon, &patchyEpsilon, sizeof(float)));
    
}

////////////////////bonded Interactions //////////////////////

__device__ void CUDAspring(c_number3 &r,c_number &r0,c_number &k,c_number4 &F){
    c_number rmod = sqrt(CUDA_DOT(r,r));
    c_number dist = abs(rmod - r0);
    F.w += 0.5f*k*dist*dist;
    c_number magForce = -(k*dist)/rmod;
    F.x -= r.x*magForce;
    F.y -= r.y*magForce;
    F.z -= r.z*magForce;
};

///////////////////nonbonded Interactions ////////////////////

///////////////////Voulume Exclusion /////////////////////////

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


///////////////////Patchy Interaction /////////////////////////

__device__ void CUDApatchy(c_number4 &r,c_number4 &pa1, c_number4 &pa2, c_number4 &pa3, c_number4 &qa1, c_number4 &qa2, c_number4 &qa3, c_number4 &F, c_number4 &tor, int ptype, int qtype){
    for(int pi = 0; pi < NumPatches[ptype]; pi++){
        c_number4 ppatch = {pa1.x * basePatchConfig[ptype][pi].x + pa2.x * basePatchConfig[ptype][pi].y + pa3.x * basePatchConfig[ptype][pi].z, pa1.y * basePatchConfig[ptype][pi].x + pa2.y * basePatchConfig[ptype][pi].y + pa3.y * basePatchConfig[ptype][pi].z, pa1.z * basePatchConfig[ptype][pi].x + pa2.z * basePatchConfig[ptype][pi].y + pa3.z * basePatchConfig[ptype][pi].z, 0};
        ppatch *= sigma;
        for(int pj = 0; pj < NumPatches[qtype]; pj++){
            c_number4 qpatch = {qa1.x * basePatchConfig[qtype][pj].x + qa2.x * basePatchConfig[qtype][pj].y + qa3.x * basePatchConfig[qtype][pj].z, qa1.y * basePatchConfig[qtype][pj].x + qa2.y * basePatchConfig[qtype][pj].y + qa3.y * basePatchConfig[qtype][pj].z, qa1.z * basePatchConfig[qtype][pj].x + qa2.z * basePatchConfig[qtype][pj].y + qa3.z * basePatchConfig[qtype][pj].z, 0};
            qpatch *= sigma;
            c_number4 patch_dist = {r.x + qpatch.x - ppatch.x, r.y + qpatch.y - ppatch.y, r.z + qpatch.z - ppatch.z, 0};
            c_number dist = CUDA_DOT(patch_dist, patch_dist);
            if(dist < patchyRcutSqr){
                c_number r8b10 = dist * dist * dist * dist / patchyAlpha;
                c_number exp_part = -1.001f * patchyEpsilon * expf(-(c_number) 0.5f * r8b10 * dist);
                F.w += exp_part;
                c_number4 tmp_force = patch_dist * (5.f * exp_part * r8b10);
                tor -= _cross(ppatch, tmp_force);
                F.x -= tmp_force.x;
                F.y -= tmp_force.y;
                F.z -= tmp_force.z;
            }
        }
    }
};

///////////////////Helix Interaction /////////////////////////

__device__ void CUDAbondedTwist(c_number4 &fp, c_number4 &fq, c_number4 &vp,c_number4 &vq,c_number4 &up, c_number4 &uq){
    c_number cos_alpha_plus_gamma = CUDA_DOT(fp,fq)+CUDA_DOT(vp,vq)/  1+CUDA_DOT(up,uq);
    
};

// __device__ void CUDAbondedDoubleBending(c_number4 &up, c_number4 &tp, c_number kb1, c_number kb2, c_number tu, c_number tk, c_number4 &F, c_number4 &tor){
//     c_number cosine = CUDA_DOT(up, uq);
//     c_number gu = cosf(tu);
// 	c_number angle = LRACOS(cosine);


//     c_number A = (kb2 * sinf(tk) - kb1 * sinf(tu)) / (tk - tu);
// 	c_number B = kb1 * sinf(tu);
// 	c_number C = kb1 * (1.f - gu) - (A * SQR(tu) / 2.f + tu * (B - tu * A));
// 	c_number D = cosf(tk) + (A * SQR(tk) * 0.5f + tk * (B - tu * A) + C) / kb2;

//     c_number g1 = (angle - tu) * A + B;
// 	c_number g_ = A * SQR(angle) / 2.f + angle * (B - tu * A);
// 	c_number g = g_ + C;

//     // unkinked bending regime
//     if(angle < tu) {
// 		tor += MD_kb[0] * kb1 * _cross(up, uq);
// 		F.w += MD_kb[0] * (1.f - cosine) * kb1;
// 	}
//     // kinked bending regime
//     else if(angle < tk) {
//         tor += MD_kb[0] * (g1 * _cross(up, uq) + g_ * _cross(up, _cross(up, uq)));
//         F.w += MD_kb[0] * (g * (1.f - cosine) + g1 * gu);
//     }
//     // unkinked bending regime
//     else {
//         tor += MD_kb[0] * kb2 * _cross(up, uq);
//         F.w += MD_kb[0] * (1.f - D) * kb2;
//     }
// }

__device__ void CUDAbondedAlignment(){
    return;
}

__device__ void CUDAbondedParticles(c_number4 &pPos, c_number4 &qPos,CUDABox *box){
    c_number4 r = box->minimum_image(pPos,qPos);
};


__device__ void CUDAnonBondedParticles(c_number4 &pPos, c_number4 &qPos, c_number4 &a1, c_number4 &a2, c_number4 &a3, c_number4 &b1, c_number4 &b2, c_number4 &b3, c_number4 &f, c_number4 &tor, CUDABox *box ){ //__device__ can only be called from a __global__ function, this only run on the GPU
    c_number4 r = box->minimum_image(pPos,qPos); //minimum distance in periodic boundary conditions
    c_number r2 = CUDA_DOT(r, r); // r^2
    if(r2 >= rcut2) return; //if r^2 is greater than total cutoff stop
}

/////////////// Main Particle Interaction //////////////////////
__global__ void CUDAparticle(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *number_neighs, CUDABox *box){
    
}




void CUDAPHBInteraction::compute_forces(CUDABaseList *lists,c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
    CUDAparticle<<<this->_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_box);
    CUT_CHECK_ERROR("Kernel failed, something quite exciting");
}
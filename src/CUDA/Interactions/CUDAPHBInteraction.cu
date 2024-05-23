#include "CUDAPHBInteraction.cuh"

__constant__ float rcut2;  // cut-off distance squared
__constant__ int exclusionType; // 0 for Linear, 1 for Cubic, 2 for Hard
__constant__ float sigma; // 
__constant__ float GPUpatchyB; // Controls the stiffness of exe volume and in case of hard the power over (sigma/r).
__constant__ float patchyRcutSqr;
__constant__ float GPUpatchyPowAlpha;
__constant__ float patchyEps;
__constant__ float hardVolCutoff;
__constant__ int GpuN;
__constant__ int n_forces;

__constant__ float patches[GPUmaxiP][5];// color,strength,x,y,z
__constant__ int numPatches[GPUmaxiC][MaxPatches];// num, patch1, patch2, patch3, patch4, patch5, patch6

__device__ int connection[MAXparticles][MAXneighbour+1];//first intger will state number of connections
__device__ float ro[MAXparticles][MAXneighbour];//rarius of the spring
__device__ float k[MAXparticles][MAXneighbour];//spring constant
__device__ int strand[MAXparticles];
__device__ int iC[MAXparticles];
__device__ float radius[MAXparticles];
__device__ float mass[MAXparticles];

CUDAPHBInteraction::CUDAPHBInteraction()
{
    _edge_compatible = true;
}

CUDAPHBInteraction::~CUDAPHBInteraction()
{
}

void CUDAPHBInteraction::get_settings(input_file &inp)
{
    PHBInteraction::get_settings(inp);
}

void CUDAPHBInteraction::cuda_init(int N)
{
    CUDABaseInteraction::cuda_init(N);
    PHBInteraction::init();

    number r8b10 = powf(patchyRcut, (number) 8.f) / patchyPowAlpha;
    number GPUhardVolCutoff = -1.001f * exp(-(number) 0.5f * r8b10 * patchyRcut2);

    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(GpuN, &N, sizeof(int)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(rcut2, &_sqr_rcut, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(sigma, &patchySigma, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(GPUpatchyB, &patchyB, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyRcutSqr, &patchyRcut2, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(GPUpatchyPowAlpha, &patchyPowAlpha, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyEps, &patchyEpsilon, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(hardVolCutoff, &GPUhardVolCutoff, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(n_forces, &_n_forces, sizeof(int)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(patches, &GPUpatches, sizeof(float) * GPUmaxiP * 5));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(numPatches, &GPUnumPatches, sizeof(int) * GPUmaxiC * MaxPatches));

    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(connection, &GPUconnections, sizeof(int) * MAXparticles * (MAXneighbour + 1));
}

////////////////////bonded Interactions //////////////////////

__device__ void CUDAspring(c_number3 &r, c_number &r0, c_number &k, c_number4 &F)
{
    c_number rmod = sqrt(CUDA_DOT(r, r));
    c_number dist = abs(rmod - r0);
    F.w += 0.5f * k * dist * dist;
    c_number magForce = -(k * dist) / rmod;
    F.x -= r.x * magForce;
    F.y -= r.y * magForce;
    F.z -= r.z * magForce;
};

///////////////////nonbonded Interactions ////////////////////

///////////////////Voulume Exclusion /////////////////////////

__device__ void CUDAexeVolLin(c_number prefactor, c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc)
{
    auto r2 = CUDA_DOT(r, r);
    if (r2 < SQR(rc))
    {
        if (r2 > SQR(rstar))
        {
            c_number rmod = sqrt(r2);
            c_number rrc = rmod - rc;
            c_number fmod = 2.f * prefactor * b * rrc / rmod;
            F.x -= r.x * fmod;
            F.y -= r.y * fmod;
            F.z -= r.z * fmod;
            F.w += prefactor * b * SQR(rrc);
            return;
        }
        c_number lj_part = CUB(SQR(sigma) / r2);
        c_number fmod = 24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / r2;
        F.x -= r.x * fmod;
        F.y -= r.y * fmod;
        F.z -= r.z * fmod;
        F.w += 4.f * prefactor * (SQR(lj_part) - lj_part);
        return;
    }
}

__device__ void CUDAexeVolCub(c_number prefactor, c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc)
{
    auto r2 = CUDA_DOT(r, r);
    if (r2 < SQR(rc))
    {
        if (r2 > SQR(rstar))
        {
            c_number rmod = sqrt(r2);
            c_number rrc = rmod - rc;
            c_number fmod = 2.f * prefactor * b * CUB(rrc) / rmod;
            F.x -= r.x * fmod;
            F.y -= r.y * fmod;
            F.z -= r.z * fmod;
            F.w += prefactor * b * SQR(SQR(rrc));
            return;
        }
        c_number lj_part = CUB(SQR(sigma) / r2);
        c_number fmod = 24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / r2;
        F.x -= r.x * fmod;
        F.y -= r.y * fmod;
        F.z -= r.z * fmod;
        F.w += 4.f * prefactor * (SQR(lj_part) - lj_part);
        return;
    }
}

__device__ void CUDAexeVolHard(c_number4 &r, c_number4 &F)
{
    c_number r2 = CUDA_DOT(r, r);
    c_number part = powf(sigma / r2, GPUpatchyB * 0.5);
    F.w += part - hardVolCutoff;
    c_number fmod = GPUpatchyB * part / r2;
    F.x -= r.x * fmod;
    F.y -= r.y * fmod;
    F.z -= r.z * fmod;
}

/////////////////// ROY Patchy Interaction /////////////////////////

// __device__ void CUDApatchy(c_number4 &r, c_number4 &pa1, c_number4 &pa2, c_number4 &pa3, c_number4 &qa1, c_number4 &qa2, c_number4 &qa3, c_number4 &F, c_number4 &tor, int ptype, int qtype)
// {
//     for (int pi = 0; pi < GPUnumPatches[ptype]; pi++)
//     {
//         c_number4 ppatch = {pa1.x * GPUbasePatchConfig[ptype][pi].x + pa2.x * GPUbasePatchConfig[ptype][pi].y + pa3.x * GPUbasePatchConfig[ptype][pi].z, pa1.y * GPUbasePatchConfig[ptype][pi].x + pa2.y * GPUbasePatchConfig[ptype][pi].y + pa3.y * GPUbasePatchConfig[ptype][pi].z, pa1.z * GPUbasePatchConfig[ptype][pi].x + pa2.z * GPUbasePatchConfig[ptype][pi].y + pa3.z * GPUbasePatchConfig[ptype][pi].z, 0};
//         ppatch *= sigma;
//         for (int pj = 0; pj < GPUnumPatches[qtype]; pj++)
//         {
//             c_number4 qpatch = {qa1.x * GPUbasePatchConfig[qtype][pj].x + qa2.x * GPUbasePatchConfig[qtype][pj].y + qa3.x * GPUbasePatchConfig[qtype][pj].z, qa1.y * GPUbasePatchConfig[qtype][pj].x + qa2.y * GPUbasePatchConfig[qtype][pj].y + qa3.y * GPUbasePatchConfig[qtype][pj].z, qa1.z * GPUbasePatchConfig[qtype][pj].x + qa2.z * GPUbasePatchConfig[qtype][pj].y + qa3.z * GPUbasePatchConfig[qtype][pj].z, 0};
//             qpatch *= sigma;
//             c_number4 patch_dist = {r.x + qpatch.x - ppatch.x, r.y + qpatch.y - ppatch.y, r.z + qpatch.z - ppatch.z, 0};
//             c_number dist = CUDA_DOT(patch_dist, patch_dist);
//             if (dist < patchyRcutSqr)
//             {
//                 c_number r8b10 = dist * dist * dist * dist / patchyAlpha;
//                 c_number exp_part = -1.001f * patchyEpsilon * expf(-0.5f * r8b10 * dist);
//                 F.w += exp_part;
//                 c_number4 tmp_force = patch_dist * (5.f * exp_part * r8b10);
//                 tor -= _cross(ppatch, tmp_force);
//                 F.x -= tmp_force.x;
//                 F.y -= tmp_force.y;
//                 F.z -= tmp_force.z;
//             }
//         }
//     }
// };


__device__ void CUDApatchyInteraction(c_number4 r, c_number4 &pa1, c_number4 &pa2, c_number4 &pa3, c_number4 &qa1, c_number4 &qa2, c_number4 &qa3, c_number4 &F, c_number4 &tor){
    
}

__device__ bool GPUbondingAllowed(int &pindex,int &qindex, int &pi, int &qi){
    // if(topology[pindex][0] == topology[qindex][0]) return false; // check if they are in the same strand
    return false;

}

///////////////////Helix Interaction /////////////////////////

__device__ void CUDAbondedTwist(c_number4 &fp, c_number4 &fq, c_number4 &vp, c_number4 &vq, c_number4 &up, c_number4 &uq)
{
    c_number cos_alpha_plus_gamma = CUDA_DOT(fp, fq) + CUDA_DOT(vp, vq) / 1 + CUDA_DOT(up, uq);
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

__device__ void CUDAbondedAlignment()
{
    return;
}

__global__ void CUDAbondedParticles(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
                                    c_number4 __restrict__ *torques, bool update_st, CUDAStressTensor *st){
};

__global__ void CUDAnonbondedParticles(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
                                       c_number4 __restrict__ *torques, const edge_bond __restrict__ *edge_list, int n_edges, bool update_st,
                                       CUDAStressTensor *st, const CUDABox *box)
{
    c_number4 dF = make_c_number4(0, 0, 0, 0);
    c_number4 dT = make_c_number4(0, 0, 0, 0);

    edge_bond b = edge_list[IND];

    // Particle 1
    c_number4 pPos = poss[b.from]; // postion
    c_number4 p1, p2, p3;          // a1,a2,a3
    get_vectors_from_quat(orientations[b.from], p1, p2, p3);

    // Particle 2
    c_number4 qPos = poss[b.to];
    c_number4 q1, q2, q3;
    get_vectors_from_quat(orientations[b.to], q1, q2, q3);

    c_number4 r = box->minimum_image(pPos, qPos);

    ////////call the main functions ////////////
    CUDAexeVolCub(2,r,dF,sigma,0.9053,GPUpatchyB,0.99998);
    // CUDApatchy(r,p1,p2,p3,q1,q2,q3,dF,dT,0,0);

    int from_index = GpuN * (IND % n_forces) + b.from;
    int to_index = GpuN * (IND % n_forces) + b.to;

    if(CUDA_DOT(dT,dT)>0.f) LR_atomicAddXYZ(&(torques[from_index]), dT);
    dT=-dT;

    if(CUDA_DOT(dF,dF)>0.f){
        LR_atomicAddXYZ(&(forces[from_index]), dF);

            if(update_st){
                CUDAStressTensor p_stress;
                _update_stress_tensor<false>(p_stress, r, dF);
                LR_atomicAddST(&(st[b.from]),p_stress);
            }

        dT += _cross(r, dF);
        LR_atomicAddXYZ(&(forces[to_index]),-dF);
    }

    if(CUDA_DOT(dT,dT)>0.f) 
    LR_atomicAddXYZ(&(torques[to_index]), dT);

}

/////////////// Main Particle Interaction //////////////////////
__global__ void CUDAparticle(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *number_neighs, CUDABox *box)
{
}

void CUDAPHBInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box)
{
    // CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);

    if(_update_st) CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
    if(_n_forces == 1){
        CUDAnonbondedParticles
        <<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
        (d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, _update_st, _d_st, d_box);
    }else{
        CUDAnonbondedParticles<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
        (d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, _update_st, _d_st, d_box);

        _sum_edge_forces_torques(d_forces, d_torques);
    }
    // CUDAbondedParticles<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
    // (d_poss, d_orientations, d_forces, d_torques, _update_st, _d_st);


    CUT_CHECK_ERROR("Kernel failed, something quite exciting");
}
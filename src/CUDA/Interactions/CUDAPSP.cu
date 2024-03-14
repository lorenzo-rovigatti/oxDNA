// Main interaction for CUDA Patchy Spring Particles


#include "CUDAPSP.cuh"

__constant__ int MD_N;
__constant__ float patches[PSPmaxPatchColor][5]; // color, strength, x, y, z // for PSP interaction patch color is useless and should be -1
__constant__ int particlePatches[PSPmaxParticleColor][PSPmaxPatchColor]; // Number of patches, patch1, patch2, patch3, patch4, patch5, patch6
// __device__ GPUconnections[PSPmaxParticles][PSPmaxNeighbour];
// __device__ float GPUr0[PSPmaxParticles][PSPmaxNeighbour]; // Equilibrium radius of the spring, if starts with 0 all particles have different radius, if 1 all particles have same radius
// __device__ float GPUk0[PSPmaxParticles][PSPmaxNeighbour]; // Spring constant, same as above
// __device__ int particleTopology[PSPmaxParticles][3]; // Strand, particleColor, radius

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

// __device__ void CUDAexeVolHard(c_number4 &r, c_number4 &F)
// {
//     c_number r2 = CUDA_DOT(r, r);
//     c_number part = powf(sigma / r2, GPUpatchyB * 0.5);
//     F.w += part - hardVolCutoff;
//     c_number fmod = GPUpatchyB * part / r2;
//     F.x -= r.x * fmod;
//     F.y -= r.y * fmod;
//     F.z -= r.z * fmod;
// }

__global__ void PSPforce(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *body3forces,
    c_number4 *torques, c_number4 *body3torques, int *matrix_neighs, int *number_neighs, cudaTextureObject_t tex_patchy_eps,
    cudaTextureObject_t tex_base_patches, bool update_st, CUDAStressTensor *st, CUDABox *box){
    
    if(IND>=MD_N) return;

    c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

    CUDAStressTensor p_st;

    }

////////////// Header Functions //////////////////////

CUDAPSPInteraction::CUDAPSPInteraction() {
};

CUDAPSPInteraction::~CUDAPSPInteraction() {
};

void CUDAPSPInteraction::get_settings(input_file &inp) {
    PSPInteraction::get_settings(inp);
}

void CUDAPSPInteraction::cuda_init(int N) {
    CUDABaseInteraction::cuda_init(N);
    PSPInteraction::init();

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));


};

void CUDAPSPInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torque, LR_bonds *d_bonds, CUDABox *d_box) {
    int N = CUDABaseInteraction::_N;
    if(_update_st) {
        CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
    }

    // PSP_forces
    //     <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
    //     (d_poss, d_orientations, d_forces, d_torque, lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _rcut, _sqr_rcut, _update_st, _d_st, d_box);
}
#include "CUDAPHBInteraction.h"

__constant__ float rcut2;  // cut-off distance squared
__constant__ int exclusionType; // 0 for Linear, 1 for Cubic, 2 for Hard
__constant__ float sigma; // 
__constant__ float patchyb; // Controls the stiffness of exe volume and in case of hard the power over (sigma/r).
__constant__ float Rstar;
__constant__ float Rc;
__constant__ float patchyRcutSqr;
__constant__ float GPUpatchyPowAlpha;
__constant__ float patchyEps;
__constant__ float hardVolCutoff;
__constant__ int MD_N[1];
__constant__ int n_forces;

__constant__ float patches[GPUmaxiP][5];// color,strength,x,y,z
__constant__ int numPatches[GPUmaxiC][MaxPatches];// num, patch1, patch2, patch3, patch4, patch5, patch6

__constant__ float ka;
__constant__ float kb;
__constant__ float kt;

__device__ int connection[MAXparticles][MAXneighbour+1];//first intger will state number of connections
__device__ float ro[MAXparticles][MAXneighbour];//rarius of the spring
__device__ float k[MAXparticles][MAXneighbour];//spring constant
__device__ int strand[MAXparticles];
__device__ int iC[MAXparticles];
__device__ float radius[MAXparticles];
__device__ float mass[MAXparticles];

// Particle properties made constant
__constant__ float tu=0.952319757;
__constant__ float tk=1.14813301;
__constant__ float kb1=1;
__constant__ float kb2=0.80;

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
    std::vector<BaseParticle *> particles(N);
	int my_N_strands;
	PHBInteraction::read_topology(&my_N_strands, particles);
    std::cout<<this->_rcut<<std::endl;
    std::cout<<this->_sqr_rcut<<std::endl;

    number r8b10 = powf(patchyRcut, (number) 8.f) / patchyPowAlpha;
    number GPUhardVolCutoff = -1.001f * exp(-(number) 0.5f * r8b10 * patchyRcut2);

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(rcut2, &_sqr_rcut, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(sigma, &patchySigma, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyb, &patchyB, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(Rstar, &patchyRstar, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(Rc, &patchyRc, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyRcutSqr, &patchyRcut2, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(GPUpatchyPowAlpha, &patchyPowAlpha, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(patchyEps, &patchyEpsilon, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(hardVolCutoff, &GPUhardVolCutoff, sizeof(float)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(n_forces, &_n_forces, sizeof(int)));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(patches, &GPUpatches, sizeof(float) * GPUmaxiP * 5));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(numPatches, &MD_NumPatches, sizeof(int) * GPUmaxiC * MaxPatches));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(connection, &GPUconnections, sizeof(int) * MAXparticles * (MAXneighbour + 1)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(ro, &GPUro, sizeof(float) * MAXparticles * MAXneighbour));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(k, &GPUk, sizeof(float) * MAXparticles * MAXneighbour));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(strand, &GPUstrand, sizeof(int) * MAXparticles));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(iC, &GPUiC, sizeof(int) * MAXparticles));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(radius, &GPUradius, sizeof(float) * MAXparticles));
    // CUDA_SAFE_CALL(cudaMemcpyToSymbol(mass, &GPUmass, sizeof(float) * MAXparticles));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(ka, &_ka, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(kb, &_kb, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(kt, &_kt, sizeof(float)));

    for(int i=0;i<7;i++){
        for(int j=0;j<GPUconnections[i][0];j++){
            std::cout<<GPUconnections[i][j+1]<<"\t";
            std::cout<<GPUro[i][j]<<"\t";
            std::cout<<GPUk[i][j]<<"\t";
        }
        std::cout<<"\n";
    }
}

__device__ void rotateVectorAroundVector(c_number4 &v, c_number4 &axis, c_number angle)
{
    // c_number4 tmp = axis * CUDA_DOT(v, axis);
    // v -= tmp;
    // c_number4 w = _cross(axis, v);
    // v = tmp + v * cosf(angle) + w * sinf(angle);

    c_number costh = cos(angle);
    c_number sinth = sin(angle);
    c_number scalar = CUDA_DOT(v, axis);
    c_number4 cross = _cross(axis, v);
    v.x = axis.x*scalar*(1.-costh)+v.x*costh+cross.x*sinth;
    v.y = axis.y*scalar*(1.-costh)+v.y*costh+cross.y*sinth;
    v.z = axis.z*scalar*(1.-costh)+v.z*costh+cross.z*sinth;
}

////////////////////bonded Interactions //////////////////////

__device__ void CUDAspring(c_number4 &r, c_number &r0, c_number &k, c_number4 &F, c_number &rmod)
{
    c_number dist = abs(rmod - r0);
    F.w -= 0.5f * k * dist * dist;
    c_number magForce = (k * dist) / rmod;
    F.x += r.x * magForce;
    F.y += r.y * magForce;
    F.z += r.z * magForce;
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
            F.x += r.x * fmod;
            F.y += r.y * fmod;
            F.z += r.z * fmod;
            F.w += prefactor * b * SQR(rrc);
        }else{
            c_number lj_part = CUB(SQR(sigma) / r2);
            c_number fmod = 24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / r2;
            F.x += r.x * fmod;
            F.y += r.y * fmod;
            F.z += r.z * fmod;
            F.w += 4.f * prefactor * (SQR(lj_part) - lj_part);
        }
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
    c_number part = powf(sigma / r2, patchyb * 0.5);
    F.w += part - hardVolCutoff;
    c_number fmod = patchyb * part / r2;
    F.x -= r.x * fmod;
    F.y -= r.y * fmod;
    F.z -= r.z * fmod;
}

/////////////////// Patchy Interaction /////////////////////////

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

// __device__ void CUDAbondedTwist(c_number4 &fp, c_number4 &fq, c_number4 &vp, c_number4 &vq, c_number4 &up, c_number4 &uq)
// {
//     c_number cos_alpha_plus_gamma = CUDA_DOT(fp, fq) + CUDA_DOT(vp, vq) / 1 + CUDA_DOT(up, uq);
// };

__device__ void CUDAbondedAlignment(c_number4 &tp, c_number4 &up,c_number4 &F, c_number4 &T,c_number &rmod){
    c_number4 force = ka*(up-tp*CUDA_DOT(up,tp)/SQR(rmod))/rmod;
    F.x += force.x;
    F.y += force.y;
    F.z += force.z;
    F.w += ka*(1-CUDA_DOT(up,tp)/rmod);
    c_number4 torque = -(ka*_cross(tp,up))/rmod;
    T.x += torque.x;
    T.y += torque.y;
    T.z += torque.z;
}

__device__ void CUDAbondedDoubleBending(c_number4 &up, c_number4 &uq,c_number4 &F, c_number4 &T){
    // c_number4 torque;

    c_number cosine = CUDA_DOT(up,uq);
    c_number gu = cosf(tu);
    c_number angle = LRACOS(cosine);

    c_number A = (kb2*sinf(tk)-kb1*sinf(tu))/(tk-tu);
    c_number B = kb1 * sinf(tu);
    c_number C = kb1 * (1.f - gu) - (A * SQR(tu) / 2.f + tu * (B - tu * A));
    c_number D = cosf(tk) + (A * SQR(tk) * 0.5f + tk * (B - tu * A) + C) / kb2;

	c_number g1 = (angle - tu) * A + B;
	c_number g_ = A * SQR(angle) / 2.f + angle * (B - tu * A);
	c_number g = g_ + C;

    if(angle < tu){
        T += kb * kb1 * _cross(up, uq);
		F.w += kb * (1.f - cosine) * kb1;
    }
    else if(angle < tk) {
		c_number4 sin_vector = _cross(up, uq);
		c_number sin_module = _module(sin_vector);
		sin_vector.x /= sin_module;
		sin_vector.y /= sin_module;
		sin_vector.z /= sin_module;

		T+= kb * g1 * sin_vector;
		F.w += kb * g;
	}
    else {
		T+= kb * kb2 * _cross(up, uq);
		F.w += kb * (D - cosine) * kb2;
	}
}

__device__ void bondedInteraction(c_number4 &poss ,c_number4 &a1, c_number4 &a2, c_number4 &a3,c_number4 &qoss, c_number4 &b1, c_number4 &b2, c_number4 &b3,c_number4 &F, c_number4 &T, CUDABox *box, int pid,int qid,c_number &ro,c_number &k){
    c_number4 r = box->minimum_image(poss, qoss);
    c_number rmod = sqrt(CUDA_DOT(r, r));
    
    CUDAspring(r,ro,k,F,rmod);
    CUDAbondedDoubleBending(a1,b1,F,T);

    if(qid-pid==1){
        r*=-1;
        CUDAbondedAlignment(r,a1,F,T,rmod);
    }else if(pid-qid==1){
        CUDAbondedAlignment(r,b1,F,T,rmod);
    }
}

__device__ void nonbondedInteraction(c_number4 &ppos ,c_number4 &a1, c_number4 &a2, c_number4 &a3,c_number4 &qpos, c_number4 &b1, c_number4 &b2, c_number4 &b3,c_number4 &F, c_number4 &T, CUDABox *box,c_number totRad){
    c_number4 r = box->minimum_image(ppos, qpos);
    
    CUDAexeVolLin(1.f,r,F,sigma*totRad,Rstar*totRad,patchyb/SQR(totRad),Rc*totRad);
    // _repulsive_lj2(1, r, F, sigma,Rstar,patchyb,Rc);
}

/////////////// Main Particle Interaction //////////////////////
__global__ void CUDAparticle(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *number_neighs, CUDABox *box){
    if(IND >= MD_N[0]) return; // for i > N leave
    c_number4 F = make_c_number4(0, 0, 0, 0);
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	// GPU_quat po = ;
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

    //Bonded Interactions
    for(int p =0;p<connection[IND][0];p++){
        int id = connection[IND][p+1];
        c_number4 b1, b2, b3;
        get_vectors_from_quat(orientations[id], b1, b2, b3);
        bondedInteraction(ppos,a1,a2,a3,poss[id],b1,b2,b3,F,T,box,IND,id,ro[IND][p],k[IND][p]);
    }


    // Non bonded interactions
    int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
    for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);
        if(k_index != IND && strand[IND]!=strand[k_index]) {
			c_number4 qpos = poss[k_index];
            c_number4 b1, b2, b3;
            get_vectors_from_quat(orientations[k_index], b1, b2, b3);
            c_number toatalRad = radius[IND] + radius[k_index];
            nonbondedInteraction(ppos,a1,a2,a3,qpos,b1,b2,b3,F,T,box,toatalRad);
        }
    }

    F.w *= (c_number) 0.5f;
    forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

void CUDAPHBInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box)
{
    // CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);

    // if(_update_st) CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
    // if(_n_forces == 1){
    //     CUDAnonbondedParticles
    //     <<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
    //     (d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, _update_st, _d_st, d_box);
    // }else{
    //     CUDAnonbondedParticles<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
    //     (d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, _update_st, _d_st, d_box);

    //     _sum_edge_forces_torques(d_forces, d_torques);
    // }
    // CUDAbondedParticles<<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
    // (d_poss, d_orientations, d_forces, d_torques, _update_st, _d_st);

    // Partiles disctributes in the grid
    CUDAparticle <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>> (d_poss, d_orientations, d_forces,  d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_box);

    
    CUT_CHECK_ERROR("Kernel failed, something quite exciting");
}
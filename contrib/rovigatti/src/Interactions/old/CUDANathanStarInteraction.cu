/*
 * CUDANathanStarInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDANathanStarInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

/* BEGIN CUDA */

/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_sqr_patchy_rcut[1];
__constant__ float MD_sqr_patchy_star_rcut[1];
__constant__ float MD_sqr_star_rcut[1];

// patch-patch constants
__constant__ int MD_rep_power[1];
__constant__ float MD_rep_E_cut[1];

__constant__ float MD_patch_angular_cutoff[1];
__constant__ int MD_patch_half_power[1];
__constant__ float MD_patch_pow_sigma[1];
__constant__ float MD_patch_pow_alpha[1];

// star-star constants
__constant__ float MD_star_f3_2[1];
__constant__ float MD_star_f1_2[1];
__constant__ float MD_T[1];
__constant__ float MD_star_sigma_s[1];
__constant__ float MD_sqr_star_sigma_s[1];
__constant__ float MD_star_factor[1];

// patchy-star variables
__constant__ int MD_interp_size[1];
__constant__ float MD_xmin[1];
__constant__ float MD_xmax[1];
__constant__ float MD_bin[1];
texture<float, 1, cudaReadModeElementType> tex_patchy_star;

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

template<typename c_number>
__device__ c_number interpolate_patchy_star(c_number x) {
	if(x <= MD_xmin[0]) return 100.;
	if(x >= MD_xmax[0]) return 0.;

	int last = (MD_interp_size[0] - 1);

	int bin_1 = (x - MD_xmin[0]) / MD_bin[0];
	int bin_0 = (bin_1 == 0) ? bin_1 : bin_1 - 1;
	int bin_2 = (bin_1 == last) ? 0 : bin_1 + 1;
	int bin_3 = (bin_2 == last) ? bin_2 : bin_2 + 1;

	c_number p0 = tex1Dfetch(tex_patchy_star, bin_0);
	c_number p1 = tex1Dfetch(tex_patchy_star, bin_1);
	c_number p2 = tex1Dfetch(tex_patchy_star, bin_2);
	c_number p3 = tex1Dfetch(tex_patchy_star, bin_3);

	c_number dx = x - bin_1 * MD_bin[0] - MD_xmin[0];
	c_number fx0 = p1;
	c_number fx1 = p2;
	c_number derx0 = (p2 - p0) / (2.f * MD_bin[0]);
	c_number derx1 = (p3 - p1) / (2.f * MD_bin[0]);
	derx0 = derx1 = 0.;

	c_number D = (2.f * (fx0 - fx1) + (derx0 + derx1) * MD_bin[0]) / MD_bin[0];
	c_number C = (fx1 - fx0 + (-derx0 - D) * MD_bin[0]);

	c_number fx = fx0 + dx * (derx0 + dx * (C + dx * D) / SQR(MD_bin[0]));

//	printf("%lf %lf %lf %lf %lf %lf\n", x, fx, p0, p1, p2, p3);
	return fx;
}

__device__ void _patchy_patchy_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &p_axis, c_number4 &q_axis, c_number4 &F, c_number4 &torque, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_patchy_rcut[0]) return;

	// here everything is done as in Allen's paper
	c_number rmod = sqrtf(sqr_r);
	c_number4 r_versor = r / (-rmod);

	// repulsion
	c_number rep_part = 1.f / powf(sqr_r, MD_rep_power[0] / 2.f);
	c_number4 force = r_versor * (MD_rep_power[0] * rep_part / rmod);
	F += force;

	c_number cospr = -CUDA_DOT(p_axis, r_versor);
	if(cospr < 0.) {
		p_axis = -p_axis;
		cospr = -cospr;
	}
	c_number cosqr = CUDA_DOT(q_axis, r_versor);
	if(cosqr < 0.) {
		q_axis = -q_axis;
		cosqr = -cosqr;
	}
	if(cospr < MD_patch_angular_cutoff[0] || cosqr < MD_patch_angular_cutoff[0]) return;

	// powf generates nan if called with a negative first argument so we have to use the SQR(...) macro and the half exponent
	c_number cospr_base = powf(SQR(cospr - 1.f), MD_patch_half_power[0] - 1) * (cospr - 1.f);
	// we do this so that later we don't have to divide this c_number by (cospr - 1), which could be 0
	c_number cospr_part = cospr_base * (cospr - 1.f);
	c_number p_mod = expf(-cospr_part / (2.f * MD_patch_pow_sigma[0]));

	c_number cosqr_base = powf(SQR(cosqr - 1.f), MD_patch_half_power[0] - 1) * (cosqr - 1.f);
	c_number cosqr_part = cosqr_base * (cosqr - 1.f);
	c_number q_mod = expf(-cosqr_part / (2.f * MD_patch_pow_sigma[0]));

	c_number sqr_surf_dist = SQR(rmod - 1.f);
	c_number r8b10 = SQR(SQR(sqr_surf_dist)) / MD_patch_pow_alpha[0];
	c_number exp_part = -1.001f * expf(-0.5f * r8b10 * sqr_surf_dist);

	// radial part
	c_number4 tmp_force = r_versor * (p_mod * q_mod * 5.f * (rmod - 1.f) * exp_part * r8b10);

	// angular p part
	c_number der_p = exp_part * q_mod * (MD_patch_half_power[0] * p_mod * cospr_base / MD_patch_pow_sigma[0]);
	c_number4 p_ortho = p_axis + cospr * r_versor;
	tmp_force -= p_ortho * (der_p / rmod);

	// angular q part
	c_number der_q = exp_part * p_mod * (-MD_patch_half_power[0] * q_mod * cosqr_base / MD_patch_pow_sigma[0]);
	c_number4 q_ortho = q_axis - cosqr * r_versor;
	tmp_force -= q_ortho * (der_q / rmod);

//	printf("%d %lf %lf %lf %lf %lf\n", IND, p_mod, q_mod, cospr_part, cosqr_part, sqrt(CUDA_DOT(tmp_force, tmp_force)));

	F += tmp_force;
	torque += _cross(r_versor, p_axis) * der_p;
}

__device__ void _patchy_star_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_patchy_star_rcut[0]) return;

	c_number mod_r = sqrt(sqr_r);
	c_number force_module = interpolate_patchy_star(mod_r);

	F.x -= r.x * (force_module / mod_r);
	F.y -= r.y * (force_module / mod_r);
	F.z -= r.z * (force_module / mod_r);
}

__device__ void _star_star_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_star_rcut[0]) return;

	c_number mod_r = sqrt(sqr_r);

	c_number common_fact = MD_star_factor[0] * 5.f * MD_T[0] * MD_star_f3_2[0] / 18.f;

	if(sqr_r < MD_sqr_star_sigma_s[0]) {
		// force over r
		c_number force_mod = common_fact / sqr_r;
		F -= r * force_mod;
	}
	else {
		c_number exp_factor = expf(-(mod_r - MD_star_sigma_s[0]) * MD_star_f1_2[0] / (2. * MD_star_sigma_s[0]));

		c_number i_f = 1.f / (1.f + MD_star_f1_2[0] * 0.5f);
		// force over r
		c_number force_mod = common_fact * i_f * exp_factor * (MD_star_sigma_s[0] / (sqr_r * mod_r) + MD_star_f1_2[0] / (2.f * sqr_r));
		F -= r * force_mod;
	}
}

// forces + second step without lists

__global__ void NS_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	int ptype = get_particle_type(ppos);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			c_number4 qpos = poss[j];
			int qtype = get_particle_type(qpos);

			int type = ptype + qtype;
			if(type == NathanStarInteraction::PATCHY_PATCHY) {
				GPU_quat qo = orientations[j];
				get_vectors_from_quat(qo, b1, b2, b3);
				_patchy_patchy_interaction(ppos, qpos, a3, b3, F, T, box);
			}
			else if(type == NathanStarInteraction::PATCHY_POLYMER) _patchy_star_interaction(ppos, qpos, F, box);
			else _star_star_interaction(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

//Forces + second step with verlet lists

__global__ void NS_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *c_number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	int ptype = get_particle_type(ppos);

	int num_neighs = c_number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		int qtype = get_particle_type(qpos);

		int type = ptype + qtype;
		if(type == NathanStarInteraction::PATCHY_PATCHY) {
			GPU_quat qo = orientations[k_index];
			get_vectors_from_quat(qo, b1, b2, b3);
			_patchy_patchy_interaction(ppos, qpos, a3, b3, F, T, box);
		}
		else if(type == NathanStarInteraction::PATCHY_POLYMER) _patchy_star_interaction(ppos, qpos, F, box);
		else _star_star_interaction(ppos, qpos, F, box);
	}

	forces[IND] = F;
	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, T);
}

/* END CUDA PART */

CUDANathanStarInteraction::CUDANathanStarInteraction() :
				CUDABaseInteraction(),
				NathanStarInteraction() {
	_d_patchy_star = NULL;
}

CUDANathanStarInteraction::~CUDANathanStarInteraction() {
	if(_d_patchy_star != NULL) CUDA_SAFE_CALL(cudaFree(_d_patchy_star));
}

void CUDANathanStarInteraction::get_settings(input_file &inp) {
	NathanStarInteraction::get_settings(inp);
}

void CUDANathanStarInteraction::_setup_cuda_interp() {
	int size = this->_spl_patchy_star->size;
	float *fx = new float[size];
	for(int i = 0; i < size; i++)
		fx[i] = -gsl_spline_eval_deriv(this->_spl_patchy_star, this->_spl_patchy_star->x[i], this->_acc_patchy_star);

	int v_size = size * sizeof(float);
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<float>(&_d_patchy_star, v_size));
	CUDA_SAFE_CALL(cudaMemcpy(_d_patchy_star, fx, size * sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_patchy_star, _d_patchy_star, v_size));

	delete[] fx;

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_interp_size, &size, sizeof(int)));
	float f_copy = this->_spl_patchy_star->interp->xmin;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_xmin, &f_copy, sizeof(float)));
	f_copy = this->_spl_patchy_star->interp->xmax;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_xmax, &f_copy, sizeof(float)));
	f_copy = this->_spl_patchy_star->x[1] - this->_spl_patchy_star->x[0];
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_bin, &f_copy, sizeof(float)));
}

void CUDANathanStarInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	NathanStarInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_rep_power, &this->_rep_power, sizeof(int)));
	int patch_half_power = this->_patch_power / 2;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_half_power, &patch_half_power, sizeof(int)));

	float f_copy = this->_sqr_patchy_rcut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_patchy_rcut, &f_copy, sizeof(float)));
	f_copy = this->_sqr_patchy_star_rcut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_patchy_star_rcut, &f_copy, sizeof(float)));
	f_copy = this->_sqr_star_rcut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_star_rcut, &f_copy, sizeof(float)));
	f_copy = this->_rep_E_cut;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_rep_E_cut, &f_copy, sizeof(float)));
	f_copy = this->_patch_angular_cutoff;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_angular_cutoff, &f_copy, sizeof(float)));
	f_copy = this->_patch_pow_alpha;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_pow_alpha, &f_copy, sizeof(float)));
	f_copy = this->_patch_pow_sigma;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_patch_pow_sigma, &f_copy, sizeof(float)));

	f_copy = this->_T;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_T, &f_copy, sizeof(float)));
	f_copy = this->_star_f1_2;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_star_f1_2, &f_copy, sizeof(float)));
	f_copy = this->_star_f3_2;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_star_f3_2, &f_copy, sizeof(float)));
	f_copy = this->_sqr_star_sigma_s;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_star_sigma_s, &f_copy, sizeof(float)));
	f_copy = this->_star_sigma_s;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_star_sigma_s, &f_copy, sizeof(float)));
	f_copy = this->_star_factor;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_star_factor, &f_copy, sizeof(float)));

	_setup_cuda_interp();
}

void CUDANathanStarInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by NathanStarInteraction");
		else {
			NS_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, d_torques, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_box);
			CUT_CHECK_ERROR("forces_second_step FS simple_lists error");
		}
	}

	CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
	if(_no_lists != NULL) {
		NS_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, d_box);
		CUT_CHECK_ERROR("forces_second_step FS no_lists error");
	}
}

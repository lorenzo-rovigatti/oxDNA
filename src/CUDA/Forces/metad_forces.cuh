/*
 * metad_forces.cuh
 *
 *  Created on: 07/may/2026
 *      Author: lorenzo
 */
#ifndef METAD_FORCES_H_
#define METAD_FORCES_H_

#include "../../Forces/Metadynamics/LTCOMTrap.h"
#include "../../Forces/Metadynamics/LTCoordination.h"

#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../../model.h"

/**
 * @brief CUDA version of a LTCOMTrap force.
 */
struct lt_com_trap {
	int type;

	int p1a_size;
	int p2a_size;
    int *p1a;
    int *p2a;

    int xmin;
    int xmax;
    int N_grid;

    c_number dX;
    c_number *potential_grid;

    int mode;
	bool PBC;
};

void init_LTCOMTrap_from_CPU(lt_com_trap *cuda_force, LTCOMTrap *cpu_force, bool first_time) {
	cuda_force->type = CUDA_LR_COM_TRAP;
	cuda_force->xmin = cpu_force->xmin;
	cuda_force->xmax = cpu_force->xmax;
	cuda_force->N_grid = cpu_force->N_grid;
	cuda_force->dX = cpu_force->dX;
	cuda_force->mode = cpu_force->_mode;
	cuda_force->PBC = cpu_force->PBC;
	cuda_force->p1a_size = cpu_force->_p1a_ptr.size();
	cuda_force->p2a_size = cpu_force->_p2a_ptr.size();

	// copy the particle arrays
	std::vector<int> local_p1a;
	for(auto particle : cpu_force->_p1a_ptr) {
		local_p1a.push_back(particle->index);
	}
	std::vector<int> local_p2a;
	for(auto particle : cpu_force->_p2a_ptr) {
		local_p2a.push_back(particle->index);
	}

	// copy the potential array
	std::vector<c_number> local_grid;
	for(auto v : cpu_force->potential_grid) {
		local_grid.push_back(v);
	}

	if(first_time) {
		// TODO: this memory never gets free'd
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&cuda_force->p1a, sizeof(int) * local_p1a.size()));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&cuda_force->p2a, sizeof(int) * local_p2a.size()));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&cuda_force->potential_grid, sizeof(c_number) * local_grid.size()));
	}
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->p1a, local_p1a.data(), sizeof(int) * local_p1a.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->p2a, local_p2a.data(), sizeof(int) * local_p2a.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->potential_grid, local_grid.data(), sizeof(c_number) * local_grid.size(), cudaMemcpyHostToDevice));
}

/**
 * @brief CUDA version of a LTCOMTrap force.
 */
struct lt_coordination {
	int type;

	int2 *pairs;
	int n_pairs;
	int idx_other;

	meta::CoordSettings::CoordMode coord_mode;
	// only used if coord_mode == MIXED, must be between 0 and 1. The contribution of the HB_ENERGY mode will be multiplied 
	// by this weight, while the contribution of the SWITCHING_FUNCTION mode will be multiplied by (1 - mixed_weight)
	c_number mixed_weight;
	c_number hb_energy_cutoff, hb_transition_width, d0, r0;
	int n;

    int coord_min;
    int coord_max;
	c_number d_coord;
    int N_grid;
    c_number *potential_grid;
};

void init_LTCoordination_from_CPU(lt_coordination *cuda_force, LTCoordination *cpu_force, int p_index, bool first_time) {
	cuda_force->type = CUDA_LR_COORDINATION;
	cuda_force->mixed_weight = cpu_force->settings.mixed_weight;
	cuda_force->coord_mode = cpu_force->settings.coord_mode;
	cuda_force->hb_energy_cutoff = cpu_force->settings.hb_energy_cutoff;
	cuda_force->hb_transition_width = cpu_force->settings.hb_transition_width;
	cuda_force->d0 = cpu_force->settings.d0;
	cuda_force->r0 = cpu_force->settings.r0;
	cuda_force->n = cpu_force->settings.n;
	cuda_force->coord_min = cpu_force->coord_min;
	cuda_force->coord_max = cpu_force->coord_max;
	cuda_force->d_coord = cpu_force->d_coord;
	cuda_force->N_grid = cpu_force->N_grid;

	// copy the pairs array
	std::vector<int2> local_pairs;
	for(auto &pair : cpu_force->all_pairs) {
		local_pairs.push_back(make_int2(pair.first->index, pair.second->index));
        if(pair.first->index == p_index) {
            cuda_force->idx_other = pair.second->index;
        }
        else if(pair.second->index == p_index) {
            cuda_force->idx_other = pair.first->index;
        }
	}
	cuda_force->n_pairs = local_pairs.size();

	// copy the potential array
	std::vector<c_number> local_grid;
	for(auto v : cpu_force->potential_grid) {
		local_grid.push_back(v);
	}

	if(first_time) {
		// TODO: this memory never gets free'd
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int2>(&cuda_force->pairs, sizeof(int2) * local_pairs.size()));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&cuda_force->potential_grid, sizeof(c_number) * local_grid.size()));
	}
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->pairs, local_pairs.data(), sizeof(int2) * local_pairs.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->potential_grid, local_grid.data(), sizeof(c_number) * local_grid.size(), cudaMemcpyHostToDevice));
}

namespace cuda_meta {

__device__ inline c_number _f1(c_number r) {
    // this should be pre-computed
	c_number shift = HYDR_EPS_OXDNA2 * SQR(1.0 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));
	
	c_number val = 0.0f;
	if(r < HYDR_RCHIGH) {
		c_number tmp = 1.0f - expf(-(r - HYDR_R0) * HYDR_A);
		val = HYDR_EPS_OXDNA2 * SQR(tmp) - shift;
	}
	return val;
}

__device__ inline c_number _f1D(c_number r) {
	c_number val = 0.0f;
	if(r < HYDR_RCHIGH) {
		c_number tmp = expf(-(r - HYDR_R0) * HYDR_A);
		val = 2.0f * HYDR_EPS_OXDNA2 * (1.0f - tmp) * tmp * HYDR_A;
	}

	return val;
}

__device__ inline c_number _f4(c_number t, c_number t0, c_number a) {
	c_number val = 0.0f;
	t -= t0;
	if(t < 0.0f) {
		t *= -1.0f;
	}
	if(t < 1.0f / sqrtf(a)) {
		val = 1.0f - a * t * t;
	}
	return val;
}

__device__ inline c_number _f4D(c_number t, c_number t0, c_number a) {
	c_number val = 0.0f;
	c_number m = 1.0f;
	c_number tt0 = t - t0;
	if(tt0 < 0.0f) {
		tt0 *= -1.0f;
		m = -1.0f;
	}
	if(tt0 < 1.0f / sqrtf(a)) {
		val = m * 2.0f * a * tt0;
	}
	return val;
}

__device__ inline c_number smooth_hb_contribution(c_number hb_energy_cutoff, c_number hb_transition_width, c_number hb_energy) {
	c_number x = (hb_energy_cutoff - hb_energy) / hb_transition_width;
	
	if(x > 10.0f) return 1.0f;
	if(x < -10.0f) return 0.0f;

	return 0.5f * (1.0f + tanhf(x));
}

__device__ inline c_number der_smooth_hb_contribution(c_number hb_energy_cutoff, c_number hb_transition_width, c_number hb_energy) {
	c_number x = (hb_energy_cutoff - hb_energy) / hb_transition_width;
	
	if(x > 10.0f) return 0.0f;
	if(x < -10.0f) return 0.0f;

	c_number tanh_x = tanhf(x);
	return -0.5f * (1.0f - SQR(tanh_x)) / hb_transition_width;
}

// Interpolate potential grid to get potential value at coordinate x
__device__ inline c_number interpolate_potential(c_number my_x, c_number dX, c_number xmin, const c_number *potential_grid, int N_grid) {
	int ix_left = (int)floorf((my_x - xmin) / dX);
	int ix_right = ix_left + 1;

	if(ix_left < 0 || ix_right >= N_grid) return 0.0f;

	c_number x_left = xmin + ix_left * dX;
	c_number w_right = (my_x - x_left) / dX;
	c_number w_left = 1.0f - w_right;

	return w_left * potential_grid[ix_left] + w_right * potential_grid[ix_right];
}

// Get derivative of potential at coordinate x
__device__ inline c_number get_x_force(c_number my_x, c_number dX, c_number xmin, const c_number *potential_grid, int N_grid) {
	int ix_left = (int)floorf((my_x - xmin) / dX);
	int ix_right = ix_left + 1;

	if(ix_left < 0 || ix_right >= N_grid) return 0.0f;

	return -(potential_grid[ix_right] - potential_grid[ix_left]) / dX;
}

// F and T are pointers so that they can be nullptr (which makes sense if compute_F_T is false)
__device__ c_number hb_interaction(const c_number4 &r, const c_number4 &rhydro, const c_number4 &ppos, const c_number4 &a1, 
    const c_number4 &a3, const c_number4 &qpos, const c_number4 &b1, const c_number4 &b3, c_number4 *F, c_number4 *T, 
    bool compute_F_T) {

    int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	c_number4 ppos_base = POS_BASE * a1;

    c_number hb_energy = 0.f;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number cost2 = -CUDA_DOT(b1, rhydrodir);
		c_number t2 = CUDA_LRACOS(cost2);
		c_number cost3 = CUDA_DOT(a1, rhydrodir);
		c_number t3 = CUDA_LRACOS(cost3);
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number cost7 = -CUDA_DOT(rhydrodir, b3);
		c_number t7 = CUDA_LRACOS(cost7);
		c_number cost8 = CUDA_DOT(rhydrodir, a3);
		c_number t8 = CUDA_LRACOS(cost8);
        
		// functions called at their relevant arguments
		c_number f1 = _f1(rhydromod);
		c_number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_A);
		c_number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_A);
		c_number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_A);
		c_number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_A);
		c_number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_A);
		c_number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_A);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if(compute_F_T && hb_energy != 0.f) {
			// derivatives called at the relevant arguments
			c_number f1D = _f1D(rhydromod);
			c_number f4t1D =  _f4D(t1, HYDR_THETA1_T0, HYDR_THETA1_A);
			c_number f4t2D =  _f4D(t2, HYDR_THETA2_T0, HYDR_THETA2_A);
			c_number f4t3D = -_f4D(t3, HYDR_THETA3_T0, HYDR_THETA3_A);
			c_number f4t4D = -_f4D(t4, HYDR_THETA4_T0, HYDR_THETA4_A);
			c_number f4t7D =  _f4D(t7, HYDR_THETA7_T0, HYDR_THETA7_A);
			c_number f4t8D = -_f4D(t8, HYDR_THETA8_T0, HYDR_THETA8_A);

			// RADIAL PART
			c_number4 Ftmp = -rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// TETA4; t4 = LRACOS (a3 * b3);
			*T += stably_normalised(_cross(a3, b3)) * (-f1 * f4t1 * f4t2 * f4t3 * f4t4D * f4t7 * f4t8);

			// TETA1; t1 = LRACOS (-a1 * b1);
			*T += stably_normalised(_cross(a1, b1)) * (-f1 * f4t1D * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
            number fact = f1 * f4t1 * f4t2D * f4t3 * f4t4 * f4t7 * f4t8;
			Ftmp += stably_normalised(b1 + rhydrodir * cost2) * (fact / rhydromod);

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
            fact = f1 * f4t1 * f4t2 * f4t3D * f4t4 * f4t7 * f4t8;
			Ftmp += stably_normalised(a1 - rhydrodir * cost3) * (fact / rhydromod);
			*T += stably_normalised(_cross(rhydrodir, a1)) * fact;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
            fact = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7D * f4t8;
			Ftmp += stably_normalised(b3 + rhydrodir * cost7) * (fact / rhydromod);

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			fact = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8D;
			Ftmp += stably_normalised(a3 - rhydrodir * cost8) * (fact / rhydromod);
			*T += stably_normalised(_cross(rhydrodir, a3)) * fact;

			*T += _cross(ppos_base, Ftmp);

			Ftmp.w = hb_energy;
			*F += Ftmp;
		}
	}

    return hb_energy;
}

__device__ c_number ipowf(c_number x, int n) {
    c_number result = 1.0f;

    while(n > 0) {
        if(n & 1) result *= x;
        x *= x;
        n >>= 1;
    }

    return result;
}

__device__ inline c_number get_pair_contribution(const lt_coordination *coord, const c_number4 &r, const c_number4 &rhydro, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a3, const c_number4 &qpos, const c_number4 &b1, const c_number4 &b3) {
	c_number coord_contribution = 0.0f;

	// Handle different coordination modes
	switch(coord->coord_mode) {
		case meta::CoordSettings::CoordMode::HB_ENERGY: {
            number hb_energy = hb_interaction(r, rhydro, ppos, a1, a3, qpos, b1, b3, nullptr, nullptr, false);
			coord_contribution = smooth_hb_contribution(coord->hb_energy_cutoff, coord->hb_transition_width, hb_energy);
			break;
		}
		case meta::CoordSettings::CoordMode::SWITCHING_FUNCTION: {
			c_number r_mod = _module(rhydro);
			c_number x = (r_mod - coord->d0) / coord->r0; // fabsf because in single precision powf can give NaNs for negative bases, even if the result should be well-defined.
			c_number xn = ipowf(x, coord->n);
			c_number switching = 1.0f / (1.0f + xn);
			coord_contribution = switching;
			break;
		}
		case meta::CoordSettings::CoordMode::MIXED: {
            number hb_energy = hb_interaction(r, rhydro, ppos, a1, a3, qpos, b1, b3, nullptr, nullptr, false);
			c_number hb_contrib = smooth_hb_contribution(coord->hb_energy_cutoff, coord->hb_transition_width, hb_energy);

			c_number r_mod = _module(rhydro);
			c_number x = (r_mod - coord->d0) / coord->r0;
			c_number xn = ipowf(x, coord->n);
			c_number switching_contrib = 1.0f / (1.0f + xn);

			coord_contribution = coord->mixed_weight * hb_contrib + (1.0f - coord->mixed_weight) * switching_contrib;
			break;
		}
	}

	return coord_contribution;
}

// Compute force and torque contribution from a single pair
// Returns force and torque contributions that need to be multiplied by df_dcoord
__device__ inline void get_pair_force_torque_contribution(const lt_coordination *coord, const c_number4 &r, 
    const c_number4 &rhydro, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a3, const c_number4 &qpos, 
    const c_number4 &b1, const c_number4 &b3, c_number4 &force, c_number4 &torque) {

    int p_idx = get_particle_index(ppos);
    int q_idx = get_particle_index(qpos);
    // Handle different coordination modes
	switch(coord->coord_mode) {
		case meta::CoordSettings::CoordMode::HB_ENERGY: {
            number hb_energy = hb_interaction(r, rhydro, ppos, a1, a3, qpos, b1, b3, &force, &torque, true);
			number dcoord_dhb_energy = der_smooth_hb_contribution(coord->hb_energy_cutoff, coord->hb_transition_width, hb_energy);
            force *= dcoord_dhb_energy;
            torque *= dcoord_dhb_energy;
			break;
		}
		case meta::CoordSettings::CoordMode::SWITCHING_FUNCTION: {
			c_number r_mod = _module(r);

			c_number x = (r_mod - coord->d0) / coord->r0;
            c_number xn = ipowf(x, coord->n);
            c_number dcoord_dr = (coord->n / coord->r0) * ipowf(x, coord->n - 1) / (SQR(1.f + xn));

            force = r / r_mod * dcoord_dr;
            torque = _cross(POS_BASE * a1, force);
			break;
		}
		case meta::CoordSettings::CoordMode::MIXED: {
            number hb_energy = hb_interaction(r, rhydro, ppos, a1, a3, qpos, b1, b3, &force, &torque, true);
			number dcoord_dhb_energy = der_smooth_hb_contribution(coord->hb_energy_cutoff, coord->hb_transition_width, hb_energy);

            c_number4 hb_force_contribution = force * dcoord_dhb_energy;
            c_number4 hb_torque_contribution = torque * dcoord_dhb_energy;

			c_number r_mod = _module(rhydro);

			c_number x = (r_mod - coord->d0) / coord->r0;
            c_number xn = ipowf(x, coord->n);
            c_number dcoord_dr = (coord->n / coord->r0) * ipowf(x, coord->n - 1) / (SQR(1.f + xn));

            c_number4 switching_force_contribution = r / r_mod * dcoord_dr;
            c_number4 switching_torque_contribution = _cross(POS_BASE * a1, switching_force_contribution);

			force = coord->mixed_weight * hb_force_contribution + (1.f - coord->mixed_weight) * switching_force_contribution;
            torque = coord->mixed_weight * hb_torque_contribution + (1.f - coord->mixed_weight) * switching_torque_contribution;
			break;
		}
	}
}

__device__ inline c_number coordination(const lt_coordination *coord, c_number4 *poss, GPU_quat *orientations, CUDABox *box) {
    number coordination = 0.f;
    for(int i = 0; i < coord->n_pairs; i++) {
        int2 pair = coord->pairs[i];
        c_number4 pos1 = poss[pair.x];
        c_number4 pos2 = poss[pair.y];

        c_number4 a1, a2, a3, b1, b2, b3;
        GPU_quat or1 = orientations[pair.x];
        GPU_quat or2 = orientations[pair.y];
        get_vectors_from_quat(or1, a1, a2, a3);
        get_vectors_from_quat(or2, b1, b2, b3);

        c_number4 r = box->minimum_image(pos1, pos2);
        c_number4 rhydro = r + POS_BASE * (b1 - a1);

        coordination += get_pair_contribution(coord, r, rhydro, pos1, a1, a3, pos2, b1, b3);
    }
    return coordination;
}

// Main function: compute force and torque for LTCoordination
// This assumes you have access to particle data (positions, orientations) in your kernel
__device__ inline void lt_coordination_force_torque(const lt_coordination *coord, int p_idx, c_number4 *poss, GPU_quat *orientations, 
    const c_number4 &ppos, c_number4 &F, c_number4 &T, CUDABox *box) {
	c_number coord_value = coordination(coord, poss, orientations, box);
    c_number df_dcoord_value = get_x_force(coord_value, coord->d_coord, coord->coord_min, coord->potential_grid, coord->N_grid);

    c_number4 qpos = poss[coord->idx_other];
    c_number4 a1, a2, a3, b1, b2, b3;
    GPU_quat por = orientations[p_idx];
    get_vectors_from_quat(por, a1, a2, a3);
    GPU_quat qor = orientations[coord->idx_other];
    get_vectors_from_quat(qor, b1, b2, b3);

    c_number4 r = box->minimum_image(ppos, qpos);
    c_number4 rhydro = r + POS_BASE * (b1 - a1);

    c_number4 force = make_c_number4(0.f, 0.f, 0.f, 0.f);
    c_number4 torque = make_c_number4(0.f, 0.f, 0.f, 0.f);
    get_pair_force_torque_contribution(coord, r, rhydro, ppos, a1, a3, qpos, b1, b3, force, torque);
    F += df_dcoord_value * force;
    T += df_dcoord_value * torque;
}

} // namespace cuda_meta

#endif /* METAD_FORCES_H_ */

#include "../cuda_utils/CUDA_lr_common.cuh"

__constant__ float MD_F1_A[2];
__constant__ float MD_F1_BLOW[2];
__constant__ float MD_F1_R0[2];
__constant__ float MD_F1_BHIGH[2];
__constant__ float MD_F1_RLOW[2];
__constant__ float MD_F1_RHIGH[2];
__constant__ float MD_F1_RCHIGH[2];
// 50 = 2 * 5 * 5
__constant__ float MD_F1_EPS[50];

__constant__ float MD_F1_RCLOW[2];
__constant__ float MD_F1_SHIFT[50];

/* System constants */
__constant__ float MD_box_side[1];

// to be included from a common file one day...
__device__ c_number4 minimum_image(const c_number4 &r_i, const c_number4 &r_j) {
	c_number dx = r_j.x - r_i.x;
	c_number dy = r_j.y - r_i.y;
	c_number dz = r_j.z - r_i.z;

	dx -= floorf(dx / MD_box_side[0] + (c_number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy / MD_box_side[0] + (c_number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz / MD_box_side[0] + (c_number) 0.5f) * MD_box_side[0];

	return make_c_number4(dx, dy, dz, (c_number) 0.f);
}

__forceinline__ __device__ c_number _f1(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		int eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BHIGH[type] * SQR(r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			c_number tmp = 1.f - expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SHIFT[eps_index];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BLOW[type] * SQR(r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _f4(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	t -= t0;
	if(t < 0) t = -t;

	if(t < tc) {
		if(t > ts) {
			// smoothing
			val = b * SQR(tc - t);
		}
		else val = (c_number) 1.f - a * SQR(t);
	}

	return val;
}

// check whether a particular pair of particles have hydrogen bonding energy lower than a given threshold hb_threshold (which may vary)
__global__ void hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb) {
	if(IND >= n_threads) return;
// for now no test for op type	if(region_is_nearhb[IND]) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = minimum_image(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets an extra two vectors that are not needed, but the function doesn't seem to be called at all, so should make little difference. 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
//	printf("\n");

	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = POS_BASE * a1;
	c_number4 qpos_base = POS_BASE * b1;

	// HYDROGEN BONDING
	c_number hb_energy = (c_number) 0;
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = _f1(rhydromod, HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	}
	// END HYDROGEN BONDING

	hb_energies[IND] = hb_energy;
}

// check whether a particular pair of particles have a 'nearly' hydrogen bond, where all or all but one of the energy factors are non-zero
__global__ void near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb) {
	if(IND >= n_threads) return;
//for now no test for op type	if(!region_is_nearhb[IND]) return; 

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = minimum_image(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets an extra two vectors that are not needed, but the function doesn't seem to be called at all, so should make little difference. 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = POS_BASE * a1;
	c_number4 qpos_base = POS_BASE * b1;

	// HYDROGEN BONDING
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);

	int total_nonzero;
	bool nearly_bonded = false;
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = _f1(rhydromod, HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		// the nearly bonded order parameter requires either all or all but one of these factors to be non-zero
		bool f1not0 = f1 < 0 ? 1 : 0;
		bool f4t1not0 = f4t1 > 0 ? 1 : 0;
		bool f4t2not0 = f4t2 > 0 ? 1 : 0;
		bool f4t3not0 = f4t3 > 0 ? 1 : 0;
		bool f4t4not0 = f4t4 > 0 ? 1 : 0;
		bool f4t7not0 = f4t7 > 0 ? 1 : 0;
		bool f4t8not0 = f4t8 > 0 ? 1 : 0;
		total_nonzero = f1not0 + f4t1not0 + f4t2not0 + f4t3not0 + f4t4not0 + f4t7not0 + f4t8not0;
		if(total_nonzero >= 6) nearly_bonded = true;
	}
	// END HYDROGEN BONDING

	nearly_bonded_array[IND] = nearly_bonded;
}

// compute the distance between a pair of particles

__global__ void dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];

	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = minimum_image(ppos, qpos);

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets an extra two vectors that are not needed, but the function doesn't seem to be called at all, so should make little difference. 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = POS_BASE * a1;
	c_number4 qpos_base = POS_BASE * b1;

	c_number4 rbase = r + qpos_base - ppos_base;
	op_dists[IND] = _module(rbase);
}

// find the minimum distances for a region and check against the stopping conditions for that region
__global__ void dist_op_eval(c_number *op_dists, int *region_lens, int *region_rows, int *cond_lens, int *cond_rows, c_number *mags, int *types, bool *stop, int n_threads, int stop_element_offset) {
	if(IND >= n_threads) return;

	// each thread checks all conditions for a given region (dimension of the order parameter)
	c_number this_value;
	int this_type;
	int this_index;
	// I'm sure this is very inefficient
	// find the minimum distance for this region
	c_number candidate;
	c_number mindist = op_dists[region_rows[IND]];
	// start at xx=1 because the xx=0 term is the starting value for mindist
	for(int xx = 1; xx < region_lens[IND]; xx++) {
		candidate = op_dists[region_rows[IND] + xx];
		if(candidate < mindist) mindist = candidate;
	}
	// check the minimum distance against each of the conditions
	for(int xx = 0; xx < cond_lens[IND]; xx++) {
		this_index = cond_rows[IND] + xx;
		stop[this_index + stop_element_offset] = false;
		this_value = mags[this_index];
		this_type = types[this_index];
		switch(this_type) {
		case 0: //GEQ
			if(mindist >= this_value) {
				stop[this_index + stop_element_offset] = true;
			}
			break;
		case 1: //LEQ
			if(mindist <= this_value) {
				stop[this_index + stop_element_offset] = true;
			}
			break;
		case 2: //LT
			if(mindist < this_value) {
				stop[this_index + stop_element_offset] = true;
			}
			break;
		case 3: //GT
			if(mindist > this_value) {
				stop[this_index + stop_element_offset] = true;
			}
			break;
		case 4: //EQ
			if(mindist == this_value) {
				stop[this_index + stop_element_offset] = true;
			}
			break;
		}
	}
}

// count up the c_number of hydrogen bonds in each region and check against the stopping conditions for that region
__global__ void hb_op_eval(float *hb_energies, int *region_lens, int *region_rows, int *cond_lens, int *cond_rows, c_number *mags, int *types, bool *stop, int n_threads, float *hb_cutoffs) {
	if(IND >= n_threads) return;

	int this_value;
	int this_type;
	int this_index;
	// find the total c_number of hbonds
	int tot_hb = 0;
	for(int xx = 0; xx < region_lens[IND]; xx++) {
		if(hb_energies[region_rows[IND] + xx] < hb_cutoffs[IND]) tot_hb += 1;
	}
	// check the c_number of hydrogen bonds against each of the conditions for this region
	for(int xx = 0; xx < cond_lens[IND]; xx++) {
		this_index = cond_rows[IND] + xx;
		stop[this_index] = false;
		this_value = mags[this_index];
		this_type = types[this_index];
		switch(this_type) {
		case 0: //GEQ
			if(tot_hb >= this_value) stop[this_index] = true;
			break;
		case 1: //LEQ
			if(tot_hb <= this_value) stop[this_index] = true;
			break;
		case 2: //LT
			if(tot_hb < this_value) stop[this_index] = true;
			break;
		case 3: //GT
			if(tot_hb > this_value) stop[this_index] = true;
			break;
		case 4: //EQ
			if(tot_hb == this_value) stop[this_index] = true;
			break;
		}
	}
}

__global__ void nearhb_op_eval(bool *nearhb_states, int *region_lens, int *region_rows, int *cond_lens, int *cond_rows, c_number *mags, int *types, bool *stop, int n_threads, int stop_element_offset) {
	if(IND >= n_threads) return;

	int this_value;
	int this_type;
	int this_index;
	// find the total c_number of hbonds
	int tot_nearhb = 0;
	for(int xx = 0; xx < region_lens[IND]; xx++) {
		if(nearhb_states[region_rows[IND] + xx]) tot_nearhb += 1;
	}
	// check the c_number of hydrogen bonds against each of the conditions for this region
	for(int xx = 0; xx < cond_lens[IND]; xx++) {
		this_index = cond_rows[IND] + xx;
		stop[this_index + stop_element_offset] = false;
		// nb implicit cast from float to int here....
		this_value = mags[this_index];
		this_type = types[this_index];
		switch(this_type) {
		case 0: //GEQ
			if(tot_nearhb >= this_value) stop[this_index + stop_element_offset] = true;
			break;
		case 1: //LEQ
			if(tot_nearhb <= this_value) stop[this_index + stop_element_offset] = true;
			break;
		case 2: //LT
			if(tot_nearhb < this_value) stop[this_index + stop_element_offset] = true;
			break;
		case 3: //GT
			if(tot_nearhb > this_value) stop[this_index + stop_element_offset] = true;
			break;
		case 4: //EQ
			if(tot_nearhb == this_value) stop[this_index + stop_element_offset] = true;
			break;
		}
	}
}


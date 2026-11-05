/*
 * CUDA_ANM.cuh
 *
 * Device code of CUDAANMInteraction: ANM/ANMT proteins together with oxDNA2, oxRNA2 and oxNA nucleotides
 * (DNANM/DNACT, RNANM/RNACT, NANM/NACT). It is the GPU version of ANMModel and ANMHybridInteraction (CPU).
 *
 * The nucleotide-nucleotide terms are those of CUDA_DRH.cuh (which must be included first), so DNA, RNA and DNA/RNA pairs are
 * computed exactly as with interaction_type = NA. The species of each particle is stored by topology index:
 * 0 for DNA, 1 for RNA and 2 for amino acids. The springs are stored by topology index too, so particles cannot be sorted
 * (CUDA_sort_every = 0).
 *
 * As on the CPU, the force of a pair term is the force acting on the second particle (the head of r = q - p).
 */

namespace cuda_anm {

enum {
	SPECIES_DNA = 0,
	SPECIES_RNA = 1,
	SPECIES_PROTEIN = 2
};

/// a quartic-smoothed repulsive Lennard-Jones (see ANMModel::RepulsiveLJ)
struct RepulsiveLJ {
	c_number sigma, rstar, b, rcut, sqr_rcut, stiffness;
};

struct ANMParameters {
	RepulsiveLJ protein, backbone, base;
	c_number protein_nucleotide_sqr_rcut;
	bool torsion;
	bool per_pair_stiffness;
	c_number k_bend, k_tor;
};

__constant__ ANMParameters anm_params;
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];

/// the length of r in the precision of the backend (_module() always uses single precision)
__forceinline__ __device__ c_number _anm_module(const c_number4 &r) {
	return sqrt(CUDA_DOT(r, r));
}

/// returns the force acting on the particle at the head of r (energy in the w component)
__forceinline__ __device__ c_number4 _anm_repulsive_lj(const c_number4 &r, const RepulsiveLJ &lj) {
	c_number4 force = make_c_number4(0, 0, 0, 0);
	c_number rnorm = CUDA_DOT(r, r);

	if(rnorm < SQR(lj.rcut)) {
		c_number fmod;
		if(rnorm > SQR(lj.rstar)) {
			c_number rmod = sqrt(rnorm);
			c_number rrc = rmod - lj.rcut;
			force.w = lj.stiffness * lj.b * SQR(SQR(rrc));
			fmod = -(4 * lj.stiffness * lj.b * CUB(rrc) / rmod);
		}
		else {
			c_number tmp = SQR(lj.sigma) / rnorm;
			c_number lj_part = tmp * tmp * tmp;
			force.w = 4 * lj.stiffness * (SQR(lj_part) - lj_part);
			fmod = -(24 * lj.stiffness * (lj_part - 2 * SQR(lj_part)) / rnorm);
		}
		force.x = r.x * fmod;
		force.y = r.y * fmod;
		force.z = r.z * fmod;
	}

	return force;
}

/**
 * @brief Bending and torsion between the consecutive amino acids p and q = p + 1 (see ANMModel::angular).
 *
 * @param r q - p
 * @param cosines equilibrium values of the four cosines of the pair (a0, b0, c0, d0)
 * @param Fq force on q (energy in the w component)
 * @param Tp, Tq torques on p and q, in the lab frame
 */
__forceinline__ __device__ void _anm_angular(const c_number4 &r, const c_number4 &a1, const c_number4 &a3, const c_number4 &b1, const c_number4 &b3,
		const c_number *cosines, c_number k_bend, c_number k_tor, c_number4 &Fq, c_number4 &Tp, c_number4 &Tq) {
	c_number rmod = _anm_module(r);
	c_number4 rij = r / rmod;
	c_number4 rji = -rij;

	c_number o1 = CUDA_DOT(rij, a1) - cosines[0];
	c_number o2 = CUDA_DOT(rji, b1) - cosines[1];
	c_number o3 = CUDA_DOT(a1, b1) - cosines[2];
	c_number o4 = CUDA_DOT(a3, b3) - cosines[3];

	Fq = -(_cross(rji, _cross(rij, a1)) * (k_bend * o1) - _cross(rji, _cross(rij, b1)) * (k_bend * o2)) / rmod;
	Fq.w = k_bend * (c_number) 0.5f * (SQR(o1) + SQR(o2)) + k_tor * (c_number) 0.5f * (SQR(o3) + SQR(o4));

	c_number4 ta = _cross(rij, a1) * (o1 * k_bend);
	c_number4 tb = _cross(rji, b1) * (o2 * k_bend);
	c_number4 torsion = _cross(a1, b1) * (o3 * k_tor) + _cross(a3, b3) * (o4 * k_tor);

	Tp = ta - torsion;
	Tq = tb + torsion;
	Tp.w = Tq.w = (c_number) 0.f;
}

/**
 * @brief Excluded volume between an amino acid and the backbone and base sites of a nucleotide (see ANMModel::nucleotide_excluded_volume).
 *
 * @param r position of the amino acid relative to the nucleotide
 * @param F_protein force on the amino acid (energy in the w component)
 * @param T_nucleotide, T_protein torques on the nucleotide and on the amino acid, in the lab frame
 */
__forceinline__ __device__ void _anm_nucleotide_excluded_volume(const c_number4 &r, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		bool nucleotide_is_rna, c_number4 &F_protein, c_number4 &T_nucleotide, c_number4 &T_protein) {
	const cuda_drh::HybridModel &hm = cuda_drh::hybrid_model;
	c_number4 back, base;
	if(nucleotide_is_rna) {
		back = a1 * hm.rna_back_a1 + a2 * hm.rna_back_a2 + a3 * hm.rna_back_a3;
		base = a1 * hm.rna_base;
	}
	else {
		back = a1 * hm.dna_back1 + a2 * hm.dna_back2;
		base = a1 * hm.dna_base;
	}

	F_protein = T_nucleotide = T_protein = make_c_number4(0, 0, 0, 0);

	const c_number4 sites[2] = { back, base };
	const RepulsiveLJ *ljs[2] = { &anm_params.backbone, &anm_params.base };
	for(int s = 0; s < 2; s++) {
		c_number4 r_site = r - sites[s];
		if(CUDA_DOT(r_site, r_site) < ljs[s]->sqr_rcut) {
			c_number4 force = _anm_repulsive_lj(r_site, *ljs[s]);
			F_protein += force;
			T_nucleotide += _cross(sites[s], -force);
			if(anm_params.torsion) {
				// point of contact on the amino acid relative to its centre
				T_protein += _cross(-stably_normalised(r_site) * (anm_params.protein.sigma * (c_number) 0.5f), force);
			}
		}
	}
	T_nucleotide.w = T_protein.w = (c_number) 0.f;
}

/**
 * @brief Nonbonded interactions for the edge-based approach. The torques are accumulated in the lab frame and moved to the frame of each
 * particle by anm_forces_edge_bonded().
 */
__global__ void anm_forces_edge_nonbonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const edge_bond __restrict__ *edge_list, int n_edges, const LR_bonds __restrict__ *bonds,
		const int * __restrict__ species, cuda_drh::KernelParameters params, const CUDABox *box) {
	if(IND >= n_edges) return;

	edge_bond b = edge_list[IND];

	c_number4 ppos = poss[b.from];
	int p_species = species[get_particle_index(ppos)];
	c_number4 qpos = poss[b.to];
	int q_species = species[get_particle_index(qpos)];

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[b.from], a1, a2, a3);
	c_number4 b1, b2, b3;
	get_vectors_from_quat(orientations[b.to], b1, b2, b3);

	c_number4 r = box->minimum_image(ppos, qpos);

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;

	bool p_protein = (p_species == SPECIES_PROTEIN);
	bool q_protein = (q_species == SPECIES_PROTEIN);

	if(!p_protein && !q_protein) {
		// same as cuda_drh::drh_forces_edge_nonbonded
		LR_bonds pbonds = bonds[b.from];
		LR_bonds qbonds = bonds[b.to];
		bool p_is_end = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);
		bool q_is_end = (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID);
		bool p_is_rna = (p_species == SPECIES_RNA);
		bool q_is_rna = (q_species == SPECIES_RNA);

		c_number4 dF = make_c_number4(0, 0, 0, 0);
		c_number4 dT = make_c_number4(0, 0, 0, 0);
		if(p_is_rna != q_is_rna) {
			cuda_drh::_particle_particle_hybrid_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, p_is_rna, p_is_end, q_is_end);
		}
		else if(p_is_rna) {
			cuda_rna::_particle_particle_RNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, params.average, true, params.mismatch_repulsion, p_is_end, q_is_end);
		}
		else {
			cuda_dna::_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, true, true, true, p_is_end, q_is_end);
		}

		c_number pair_energy = dF.w;
		c_number hb_energy = dT.w;
		dF.w = (c_number) 0;
		dT.w = (c_number) 0;

		if(CUDA_DOT(dT, dT) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[from_index]), dT);
		dT = -dT; // torque acting on particle q (which also has another contribution if dF > 0)

		if(CUDA_DOT(dF, dF) > (c_number) 0.f) {
			LR_atomicAddXYZ(&(forces[from_index]), dF);
			// Allen Eq. 6 pag 3:
			dT += _cross(r, dF);
			LR_atomicAddXYZ(&(forces[to_index]), -dF);
		}

		if(CUDA_DOT(dT, dT) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[to_index]), dT);

		if(pair_energy != (c_number) 0) {
			atomicAdd(&(forces[from_index].w), pair_energy);
			atomicAdd(&(forces[to_index].w), pair_energy);
		}
		if(hb_energy != (c_number) 0) {
			atomicAdd(&(torques[from_index].w), hb_energy);
			atomicAdd(&(torques[to_index].w), hb_energy);
		}
		return;
	}

	if(p_protein && q_protein) {
		// pairs connected by n3/n5 (backbone neighbours with ANMT) are not in the edge list: see anm_forces_edge_bonded
		if(CUDA_DOT(r, r) >= anm_params.protein.sqr_rcut) return;
		c_number4 force = _anm_repulsive_lj(r, anm_params.protein);
		c_number energy = force.w;
		force.w = (c_number) 0.f;
		if(CUDA_DOT(force, force) > (c_number) 0.f) {
			LR_atomicAddXYZ(&(forces[from_index]), -force);
			LR_atomicAddXYZ(&(forces[to_index]), force);
			if(anm_params.torsion) {
				c_number4 half_r = r * (c_number) 0.5f;
				LR_atomicAddXYZ(&(torques[from_index]), _cross(half_r, -force));
				LR_atomicAddXYZ(&(torques[to_index]), _cross(-half_r, force));
			}
		}
		if(energy != (c_number) 0) {
			atomicAdd(&(forces[from_index].w), energy / 2);
			atomicAdd(&(forces[to_index].w), energy / 2);
		}
		return;
	}

	// amino acid - nucleotide
	if(CUDA_DOT(r, r) >= anm_params.protein_nucleotide_sqr_rcut) return;

	int protein_index, nucleotide_index;
	c_number4 r_pn;
	c_number4 n1, n2, n3;
	bool nucleotide_is_rna;
	if(p_protein) {
		protein_index = from_index;
		nucleotide_index = to_index;
		r_pn = -r;
		n1 = b1;
		n2 = b2;
		n3 = b3;
		nucleotide_is_rna = (q_species == SPECIES_RNA);
	}
	else {
		protein_index = to_index;
		nucleotide_index = from_index;
		r_pn = r;
		n1 = a1;
		n2 = a2;
		n3 = a3;
		nucleotide_is_rna = (p_species == SPECIES_RNA);
	}

	c_number4 F_protein, T_nucleotide, T_protein;
	_anm_nucleotide_excluded_volume(r_pn, n1, n2, n3, nucleotide_is_rna, F_protein, T_nucleotide, T_protein);
	c_number energy = F_protein.w;
	F_protein.w = (c_number) 0.f;

	if(CUDA_DOT(F_protein, F_protein) > (c_number) 0.f) {
		LR_atomicAddXYZ(&(forces[protein_index]), F_protein);
		LR_atomicAddXYZ(&(forces[nucleotide_index]), -F_protein);
	}
	if(CUDA_DOT(T_nucleotide, T_nucleotide) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[nucleotide_index]), T_nucleotide);
	if(CUDA_DOT(T_protein, T_protein) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[protein_index]), T_protein);
	if(energy != (c_number) 0) {
		atomicAdd(&(forces[protein_index].w), energy / 2);
		atomicAdd(&(forces[nucleotide_index].w), energy / 2);
	}
}

/**
 * @brief Bonded interactions of an amino acid: springs and, with ANMT, the excluded volume, bending and torsion with its backbone neighbours
 * (n3 = previous, n5 = next amino acid), which are not in the neighbour lists. Adds the force on p to F and the torque (lab frame) to T.
 *
 * The springs of the particle with topology index i are those in [spring_start[i], spring_start[i + 1]) of spring_neighbour,
 * spring_k and spring_r0. The bending/torsion parameters of the pair (i, i + 1) are the four cosines angular_cosines[4 * i] ... [4 * i + 3] and, if they are
 * given per pair, the stiffnesses angular_k[2 * i] (bending) and angular_k[2 * i + 1] (torsion).
 */
__forceinline__ __device__ void _anm_protein_bonded(const c_number4 *poss, const GPU_quat *orientations, const LR_bonds &pbonds, const c_number4 &ppos,
		int p_index, const c_number4 &a1, const c_number4 &a3, const int *spring_start, const int *spring_neighbour, const c_number *spring_k,
		const c_number *spring_r0, const c_number *angular_cosines, const c_number *angular_k, const CUDABox *box, c_number4 &F, c_number4 &T) {
	// springs. Particles are not sorted, so topology indices are also the indices of the arrays
	for(int i = spring_start[p_index]; i < spring_start[p_index + 1]; i++) {
		c_number4 r = box->minimum_image(ppos, poss[spring_neighbour[i]]);
		c_number d = _anm_module(r);
		c_number disp = d - spring_r0[i];
		// force on p (the CPU code computes the force on q, -k * disp / d * r)
		c_number4 dF = r * (spring_k[i] * disp / d);
		dF.w = (c_number) 0.25f * spring_k[i] * SQR(disp); // half of the energy of the pair
		F += dF;
	}

	if(anm_params.torsion) {
		// the backbone neighbours are not in the neighbour lists, so their excluded volume is computed here
		int neighbours[2] = { pbonds.n3, pbonds.n5 };
		for(int n = 0; n < 2; n++) {
			if(neighbours[n] == P_INVALID) continue;
			c_number4 r = box->minimum_image(ppos, poss[neighbours[n]]);
			c_number4 force = _anm_repulsive_lj(r, anm_params.protein);
			c_number energy = force.w;
			force.w = (c_number) 0.f;
			F -= force;
			F.w += energy / 2;
			T += _cross(r * (c_number) 0.5f, -force);
		}

		// bending and torsion with the next amino acid: p is the first particle of the pair
		if(pbonds.n5 != P_INVALID && get_particle_index(poss[pbonds.n5]) == p_index + 1) {
			c_number4 qpos = poss[pbonds.n5];
			c_number4 b1, b2, b3;
			get_vectors_from_quat(orientations[pbonds.n5], b1, b2, b3);
			c_number k_bend = anm_params.k_bend, k_tor = anm_params.k_tor;
			if(anm_params.per_pair_stiffness) {
				k_bend = angular_k[2 * p_index];
				k_tor = angular_k[2 * p_index + 1];
			}
			c_number4 Fq, Tp, Tq;
			_anm_angular(box->minimum_image(ppos, qpos), a1, a3, b1, b3, angular_cosines + 4 * p_index, k_bend, k_tor, Fq, Tp, Tq);
			c_number energy = Fq.w;
			Fq.w = (c_number) 0.f;
			F -= Fq;
			F.w += energy / 2;
			T += Tp;
		}

		// bending and torsion with the previous amino acid: p is the second particle of the pair
		if(pbonds.n3 != P_INVALID && get_particle_index(poss[pbonds.n3]) == p_index - 1) {
			c_number4 opos = poss[pbonds.n3];
			c_number4 b1, b2, b3;
			get_vectors_from_quat(orientations[pbonds.n3], b1, b2, b3);
			int o_index = p_index - 1;
			c_number k_bend = anm_params.k_bend, k_tor = anm_params.k_tor;
			if(anm_params.per_pair_stiffness) {
				k_bend = angular_k[2 * o_index];
				k_tor = angular_k[2 * o_index + 1];
			}
			c_number4 Fq, To, Tq;
			_anm_angular(box->minimum_image(opos, ppos), b1, b3, a1, a3, angular_cosines + 4 * o_index, k_bend, k_tor, Fq, To, Tq);
			c_number energy = Fq.w;
			Fq.w = (c_number) 0.f;
			F += Fq;
			F.w += energy / 2;
			T += Tq;
		}
	}
}

/**
 * @brief Bonded interactions for the edge-based approach: nucleotides get the bonded terms of their nucleic acid model, amino acids those of
 * _anm_protein_bonded(). It also changes the frame of the torques from the one of the lab to the one of each particle.
 */
__global__ void anm_forces_edge_bonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const LR_bonds __restrict__ *bonds, const int * __restrict__ species, cuda_drh::KernelParameters params,
		const int * __restrict__ spring_start, const int * __restrict__ spring_neighbour, const c_number * __restrict__ spring_k,
		const c_number * __restrict__ spring_r0, const c_number * __restrict__ angular_cosines, const c_number * __restrict__ angular_k,
		const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];
	int p_index = get_particle_index(ppos);
	int p_species = species[p_index];

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(p_species != SPECIES_PROTEIN) {
		cuda_drh::_drh_bonded_part(poss, orientations, pbonds, ppos, a1, a2, a3, p_species == SPECIES_RNA, params, F, T);
	}
	else {
		_anm_protein_bonded(poss, orientations, pbonds, ppos, p_index, a1, a3, spring_start, spring_neighbour, spring_k, spring_r0, angular_cosines,
				angular_k, box, F, T);
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

/**
 * @brief Forces and torques with neighbour matrices (use_edge = false): each thread computes the force and torque acting on its own particle.
 */
__global__ void anm_forces(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs, const int *number_neighs, const LR_bonds __restrict__ *bonds,
		const int * __restrict__ species, cuda_drh::KernelParameters params, const int * __restrict__ spring_start,
		const int * __restrict__ spring_neighbour, const c_number * __restrict__ spring_k, const c_number * __restrict__ spring_r0,
		const c_number * __restrict__ angular_cosines, const c_number * __restrict__ angular_k, const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];
	int p_index = get_particle_index(ppos);
	int p_species = species[p_index];
	bool p_protein = (p_species == SPECIES_PROTEIN);
	bool p_is_rna = (p_species == SPECIES_RNA);
	bool p_is_end = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(p_protein) {
		_anm_protein_bonded(poss, orientations, pbonds, ppos, p_index, a1, a3, spring_start, spring_neighbour, spring_k, spring_r0, angular_cosines,
				angular_k, box, F, T);
	}
	else {
		cuda_drh::_drh_bonded_part(poss, orientations, pbonds, ppos, a1, a2, a3, p_is_rna, params, F, T);
	}

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	T.w = (c_number) 0;
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);
		if(k_index == IND || pbonds.n3 == k_index || pbonds.n5 == k_index) continue;

		c_number4 qpos = poss[k_index];
		int q_species = species[get_particle_index(qpos)];
		bool q_protein = (q_species == SPECIES_PROTEIN);
		c_number4 r = box->minimum_image(ppos, qpos);
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[k_index], b1, b2, b3);

		if(!p_protein && !q_protein) {
			// same as cuda_drh::drh_forces
			LR_bonds qbonds = bonds[k_index];
			bool q_is_end = (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID);
			bool q_is_rna = (q_species == SPECIES_RNA);
			if(p_is_rna != q_is_rna) {
				cuda_drh::_particle_particle_hybrid_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, p_is_rna, p_is_end, q_is_end);
			}
			else if(p_is_rna) {
				cuda_rna::_particle_particle_RNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, params.average, true, params.mismatch_repulsion, p_is_end, q_is_end);
			}
			else {
				c_number4 dF = make_c_number4(0, 0, 0, 0);
				cuda_dna::_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, T, true, true, true, p_is_end, q_is_end);
				F += dF;
			}
		}
		else if(p_protein && q_protein) {
			if(CUDA_DOT(r, r) >= anm_params.protein.sqr_rcut) continue;
			// force on q
			c_number4 force = _anm_repulsive_lj(r, anm_params.protein);
			c_number energy = force.w;
			force.w = (c_number) 0.f;
			F -= force;
			F.w += energy / 2;
			if(anm_params.torsion) {
				T += _cross(r * (c_number) 0.5f, -force);
			}
		}
		else {
			if(CUDA_DOT(r, r) >= anm_params.protein_nucleotide_sqr_rcut) continue;
			c_number4 F_protein, T_nucleotide, T_protein;
			if(p_protein) {
				// r_pn = p - q
				_anm_nucleotide_excluded_volume(-r, b1, b2, b3, q_species == SPECIES_RNA, F_protein, T_nucleotide, T_protein);
				c_number energy = F_protein.w;
				F_protein.w = (c_number) 0.f;
				F += F_protein;
				F.w += energy / 2;
				T += T_protein;
			}
			else {
				_anm_nucleotide_excluded_volume(r, a1, a2, a3, p_is_rna, F_protein, T_nucleotide, T_protein);
				c_number energy = F_protein.w;
				F_protein.w = (c_number) 0.f;
				F -= F_protein;
				F.w += energy / 2;
				T += T_nucleotide;
			}
		}
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

} // namespace cuda_anm

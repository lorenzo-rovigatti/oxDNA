/*
 * CUDA_DRH.cuh
 *
 * Device code of CUDADRHInteraction (interaction_type = NA, i.e. oxDNA2 + oxRNA2 + DNA/RNA hybrid interactions).
 *
 * The DNA-DNA and RNA-RNA terms are the ones of CUDA_DNA.cuh (oxDNA2) and CUDA_RNA.cuh (oxRNA2), which are included in
 * their own namespaces so that their (identically named) constants and helper functions do not clash. Only the terms that
 * act between a DNA and an RNA nucleotide (excluded volume, hydrogen bonding, cross stacking, coaxial stacking and
 * Debye-Huckel), which in the CPU code are implemented by DRHInteraction, are new and live in namespace cuda_drh.
 *
 * This file is included by more than one translation unit (CUDADRHInteraction.cu and CUDAANMInteraction.cu), so everything
 * is in an anonymous namespace: each of them gets its own copy of the kernels and of the constants, which it initialises
 * itself (see CUDA_NA_constants.cuh). Without it the host stubs of the kernels are defined twice and the library does not link.
 */

#ifndef CUDA_DRH_CUH_
#define CUDA_DRH_CUH_

#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../../Interactions/drh_model.h"

namespace {

namespace cuda_rna {
#include "CUDA_RNA.cuh"
}

namespace cuda_dna {
#include "CUDA_DNA.cuh"
}

namespace cuda_drh {

/// number of nucleotides
__constant__ int MD_N[1];
/// number of force slots per nucleotide (edge_n_forces). Used only if use_edge = true
__constant__ int MD_n_forces[1];

/**
 * @brief Parameters of the DNA/RNA hybrid interaction. They are stored with the precision of the CUDA backend (c_number) and not as floats as
 * in CUDA_DNA.cuh and CUDA_RNA.cuh, so that the new terms are as precise as the backend allows.
 */
struct HybridModel {
	// positions of the interaction sites in the reference frame of a nucleotide, for a DNA (oxDNA2, with grooving) and for an RNA nucleotide
	c_number dna_back1, dna_back2, dna_base, dna_stack;
	c_number rna_back_a1, rna_back_a2, rna_back_a3, rna_base, rna_stack;

	// radial part of the hydrogen bonding: the eps and shift tables are indexed by n3 * 5 + n5
	c_number f1_eps[25], f1_shift[25];
	c_number f1_a, f1_r0, f1_blow, f1_bhigh, f1_rlow, f1_rhigh, f1_rclow, f1_rchigh;
	c_number hb_multiplier;

	// radial part of the cross (type 0) and coaxial (type 1) stacking
	c_number f2_k[2], f2_r0[2], f2_rc[2], f2_blow[2], f2_bhigh[2], f2_rlow[2], f2_rhigh[2], f2_rclow[2], f2_rchigh[2];

	// Debye-Huckel: the hybrid pairs use the oxRNA2 parameters
	c_number dh_rc, dh_rhigh, dh_b, dh_minus_kappa, dh_prefactor;
	bool dh_half_charged_ends;
};

__constant__ HybridModel hybrid_model;

/// parameters of an f4 modulation: f4(t) = 1 - a * (t - t0)^2 for |t - t0| < ts, b * (tc - |t - t0|)^2 for ts < |t - t0| < tc, 0 otherwise
struct F4Parameters {
	c_number t0, ts, tc, a, b;
};

#define DRH_F4_PARAMETERS(NAME) { (c_number) NAME##_T0, (c_number) NAME##_TS, (c_number) NAME##_TC, (c_number) NAME##_A, (c_number) NAME##_B }

// The hybrid hydrogen-bonding and cross-stacking angular modulations are the (mesh-interpolated) ones of the DNA model, which are
// numerically identical to the DRH_* constants of drh_model.h
__constant__ F4Parameters HYDR_THETA1 = DRH_F4_PARAMETERS(DRH_HYDR_THETA1);
__constant__ F4Parameters HYDR_THETA2 = DRH_F4_PARAMETERS(DRH_HYDR_THETA2);
__constant__ F4Parameters HYDR_THETA4 = DRH_F4_PARAMETERS(DRH_HYDR_THETA4);
__constant__ F4Parameters HYDR_THETA7 = DRH_F4_PARAMETERS(DRH_HYDR_THETA7);
__constant__ F4Parameters CRST_THETA1 = DRH_F4_PARAMETERS(DRH_CRST_THETA1);
__constant__ F4Parameters CRST_THETA2 = DRH_F4_PARAMETERS(DRH_CRST_THETA2);
__constant__ F4Parameters CRST_THETA4 = DRH_F4_PARAMETERS(DRH_CRST_THETA4);
__constant__ F4Parameters CRST_THETA7 = DRH_F4_PARAMETERS(DRH_CRST_THETA7);
// as in DRHInteraction, the hybrid coaxial stacking uses the oxDNA (not oxDNA2) value of the centre of f4(theta1)
__constant__ F4Parameters CXST_THETA1 = { (c_number) DRH_CXST_THETA1_T0_OXDNA, (c_number) DRH_CXST_THETA1_TS, (c_number) DRH_CXST_THETA1_TC, (c_number) DRH_CXST_THETA1_A, (c_number) DRH_CXST_THETA1_B };
__constant__ F4Parameters CXST_THETA4 = DRH_F4_PARAMETERS(DRH_CXST_THETA4);
__constant__ F4Parameters CXST_THETA5 = DRH_F4_PARAMETERS(DRH_CXST_THETA5);

__forceinline__ __device__ c_number _f4(c_number t, const F4Parameters &p) {
	c_number val = 0;
	t = fabs(t - p.t0);
	if(t < p.tc) {
		val = (t > p.ts) ? p.b * SQR(p.tc - t) : 1 - p.a * SQR(t);
	}
	return val;
}

/// derivative of f4 with respect to its argument divided by sin(t)
__forceinline__ __device__ c_number _f4Dsin(c_number t, const F4Parameters &p) {
	c_number val = 0;
	c_number tt0 = t - p.t0;
	// f4 is a parabola centred in t0: on the other side of t0 its value is the same but its derivative has the opposite sign
	c_number m = (tt0 < 0) ? -1 : 1;
	tt0 = fabs(tt0);
	if(tt0 < p.tc) {
		c_number sint = sin(t);
		if(tt0 > p.ts) {
			val = m * 2 * p.b * (tt0 - p.tc) / sint;
		}
		else {
			val = (SQR(sint) > (c_number) 1e-8) ? -m * 2 * p.a * tt0 / sint : -m * 2 * p.a;
		}
	}
	return val;
}

/// f4 and its derivative with respect to cos(t), as a function of cos(t)
__forceinline__ __device__ void _f4_of_cos(c_number cost, const F4Parameters &p, c_number &f4, c_number &df4_dcos) {
	c_number t = CUDA_LRACOS(cost);
	f4 = _f4(t, p);
	df4_dcos = -_f4Dsin(t, p);
}

/// the sum f4(cos(t)) + f4(-cos(t)) and its derivative with respect to cos(t)
__forceinline__ __device__ void _f4_of_cos_symmetrised(c_number cost, const F4Parameters &p, c_number &f4, c_number &df4_dcos) {
	c_number f4_plus, d_plus, f4_minus, d_minus;
	_f4_of_cos(cost, p, f4_plus, d_plus);
	_f4_of_cos(-cost, p, f4_minus, d_minus);
	f4 = f4_plus + f4_minus;
	df4_dcos = d_plus - d_minus;
}

__forceinline__ __device__ c_number _f5(c_number f, c_number a, c_number b, c_number xc, c_number xs) {
	c_number val = 0;
	if(f > xc) {
		if(f < xs) {
			val = b * SQR(xc - f);
		}
		else if(f < 0) {
			val = 1 - a * SQR(f);
		}
		else {
			val = 1;
		}
	}
	return val;
}

__forceinline__ __device__ c_number _f5D(c_number f, c_number a, c_number b, c_number xc, c_number xs) {
	c_number val = 0;
	if(f > xc) {
		if(f < xs) {
			val = 2 * b * (f - xc);
		}
		else if(f < 0) {
			val = -2 * a * f;
		}
	}
	return val;
}

/// maps the types of the two nucleotides to the entry of the hybrid hydrogen-bonding table: AU pairs use the AA slot, and the strength of GC
/// pairs depends on which of the two nucleotides is the RNA one
__forceinline__ __device__ int _hybrid_f1_index(int n3, int n5, bool is_n3_DNA, bool is_n5_DNA) {
	if((!is_n3_DNA && n3 == N_T && n5 == N_A) || (!is_n5_DNA && n5 == N_T && n3 == N_A)) {
		n3 = n5 = N_A;
	}

	if((n3 == N_G && n5 == N_C) || (n5 == N_G && n3 == N_C)) {
		if((!is_n3_DNA && n3 == N_G) || (!is_n5_DNA && n5 == N_G)) {
			n3 = N_G;
			n5 = N_C;
		}
		else {
			n3 = N_C;
			n5 = N_G;
		}
	}

	return n3 * 5 + n5;
}

__forceinline__ __device__ c_number _hybrid_f1(c_number r, int index) {
	const HybridModel &m = hybrid_model;
	c_number val = 0;
	if(r < m.f1_rchigh) {
		if(r > m.f1_rhigh) {
			val = m.f1_eps[index] * m.f1_bhigh * SQR(r - m.f1_rchigh);
		}
		else if(r > m.f1_rlow) {
			c_number tmp = 1 - exp(-(r - m.f1_r0) * m.f1_a);
			val = m.f1_eps[index] * SQR(tmp) - m.f1_shift[index];
		}
		else if(r > m.f1_rclow) {
			val = m.f1_eps[index] * m.f1_blow * SQR(r - m.f1_rclow);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _hybrid_f1D(c_number r, int index) {
	const HybridModel &m = hybrid_model;
	c_number val = 0;
	if(r < m.f1_rchigh) {
		if(r > m.f1_rhigh) {
			val = 2 * m.f1_bhigh * (r - m.f1_rchigh);
		}
		else if(r > m.f1_rlow) {
			c_number tmp = exp(-(r - m.f1_r0) * m.f1_a);
			val = 2 * (1 - tmp) * tmp * m.f1_a;
		}
		else if(r > m.f1_rclow) {
			val = 2 * m.f1_blow * (r - m.f1_rclow);
		}
	}
	return m.f1_eps[index] * val;
}

__forceinline__ __device__ c_number _hybrid_f2(c_number r, int type) {
	const HybridModel &m = hybrid_model;
	c_number val = 0;
	if(r < m.f2_rchigh[type]) {
		if(r > m.f2_rhigh[type]) {
			val = m.f2_k[type] * m.f2_bhigh[type] * SQR(r - m.f2_rchigh[type]);
		}
		else if(r > m.f2_rlow[type]) {
			val = (m.f2_k[type] / 2) * (SQR(r - m.f2_r0[type]) - SQR(m.f2_rc[type] - m.f2_r0[type]));
		}
		else if(r > m.f2_rclow[type]) {
			val = m.f2_k[type] * m.f2_blow[type] * SQR(r - m.f2_rclow[type]);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _hybrid_f2D(c_number r, int type) {
	const HybridModel &m = hybrid_model;
	c_number val = 0;
	if(r < m.f2_rchigh[type]) {
		if(r > m.f2_rhigh[type]) {
			val = 2 * m.f2_k[type] * m.f2_bhigh[type] * (r - m.f2_rchigh[type]);
		}
		else if(r > m.f2_rlow[type]) {
			val = m.f2_k[type] * (r - m.f2_r0[type]);
		}
		else if(r > m.f2_rclow[type]) {
			val = 2 * m.f2_k[type] * m.f2_blow[type] * (r - m.f2_rclow[type]);
		}
	}
	return val;
}

/// repulsive Lennard-Jones, smoothed beyond rstar. On return F is the force acting on the first particle and F.w the energy
__forceinline__ __device__ void _hybrid_excluded_volume(const c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc) {
	c_number rsqr = CUDA_DOT(r, r);

	F.x = F.y = F.z = F.w = 0;
	if(rsqr < SQR(rc)) {
		c_number fmod, energy;
		if(rsqr > SQR(rstar)) {
			c_number rmod = sqrt(rsqr);
			c_number rrc = rmod - rc;
			fmod = 2 * DRH_EXCL_EPS * b * rrc / rmod;
			energy = DRH_EXCL_EPS * b * SQR(rrc);
		}
		else {
			c_number lj_part = CUB(SQR(sigma) / rsqr);
			fmod = 24 * DRH_EXCL_EPS * (lj_part - 2 * SQR(lj_part)) / rsqr;
			energy = 4 * DRH_EXCL_EPS * (SQR(lj_part) - lj_part);
		}
		F.x = r.x * fmod;
		F.y = r.y * fmod;
		F.z = r.z * fmod;
		F.w = energy;
	}
}

/// positions of the backbone, base and stacking sites of a DNA or RNA nucleotide, in the lab frame
__forceinline__ __device__ void _hybrid_sites(bool is_rna, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3, c_number4 &back, c_number4 &base, c_number4 &stack) {
	const HybridModel &m = hybrid_model;
	if(is_rna) {
		back = a1 * m.rna_back_a1 + a2 * m.rna_back_a2 + a3 * m.rna_back_a3;
		base = a1 * m.rna_base;
		stack = a1 * m.rna_stack;
	}
	else {
		back = a1 * m.dna_back1 + a2 * m.dna_back2;
		base = a1 * m.dna_base;
		stack = a1 * m.dna_stack;
	}
}

/**
 * @brief Computes the interaction between nucleotide p and nucleotide q, one of which is DNA and the other RNA.
 *
 * It is the transcription of DRHInteraction::_nonbonded_excluded_volume_DRH, _hydrogen_bonding_DRH, _cross_stacking_DRH,
 * _coaxial_stacking_DRH and of the oxRNA2 Debye-Huckel term, from the point of view of p. Like the DNA and RNA functions, it adds
 * the force on p (and the energy in F.w) to F, the torque on p in the lab frame to T and the hydrogen-bonding energy to T.w.
 *
 * @param r minimum-image vector from p to q
 * @param p_is_rna true if p is the RNA nucleotide, false if p is the DNA one (q is the other one)
 */
__device__ void _particle_particle_hybrid_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1, const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, bool p_is_rna, bool p_is_end,
		bool q_is_end) {
	const HybridModel &m = hybrid_model;
	const int ptype = get_particle_type(ppos);
	const int qtype = get_particle_type(qpos);
	const int pbtype = get_particle_btype(ppos);
	const int qbtype = get_particle_btype(qpos);

	c_number4 ppos_back, ppos_base, ppos_stack, qpos_back, qpos_base, qpos_stack;
	_hybrid_sites(p_is_rna, a1, a2, a3, ppos_back, ppos_base, ppos_stack);
	_hybrid_sites(!p_is_rna, b1, b2, b3, qpos_back, qpos_base, qpos_stack);

	c_number tot_Tw = T.w;
	c_number4 Ftmp, Ttmp;

	// EXCLUDED VOLUME
	// BACK-BACK
	c_number4 rbackbone = r + qpos_back - ppos_back;
	_hybrid_excluded_volume(rbackbone, Ftmp, DRH_EXCL_S1, DRH_EXCL_R1, DRH_EXCL_B1, DRH_EXCL_RC1);
	Ttmp = _cross(ppos_back, Ftmp);
	F += Ftmp;
	// BASE-BASE
	c_number4 rhydro = r + qpos_base - ppos_base;
	_hybrid_excluded_volume(rhydro, Ftmp, DRH_EXCL_S2, DRH_EXCL_R2, DRH_EXCL_B2, DRH_EXCL_RC2);
	Ttmp += _cross(ppos_base, Ftmp);
	F += Ftmp;
	// p-BACK vs q-BASE
	_hybrid_excluded_volume(r + qpos_base - ppos_back, Ftmp, DRH_EXCL_S4, DRH_EXCL_R4, DRH_EXCL_B4, DRH_EXCL_RC4);
	Ttmp += _cross(ppos_back, Ftmp);
	F += Ftmp;
	// p-BASE vs q-BACK
	_hybrid_excluded_volume(r + qpos_back - ppos_base, Ftmp, DRH_EXCL_S3, DRH_EXCL_R3, DRH_EXCL_B3, DRH_EXCL_RC3);
	Ttmp += _cross(ppos_base, Ftmp);
	F += Ftmp;

	// DEBYE-HUCKEL
	c_number rbackmod = _module(rbackbone);
	if(rbackmod < m.dh_rc) {
		c_number4 rbackdir = rbackbone / rbackmod;
		if(rbackmod < m.dh_rhigh) {
			c_number exp_part = exp(m.dh_minus_kappa * rbackmod);
			Ftmp = rbackdir * (-m.dh_prefactor * exp_part * (m.dh_minus_kappa / rbackmod - 1 / SQR(rbackmod)));
			Ftmp.w = exp_part * (m.dh_prefactor / rbackmod);
		}
		else {
			Ftmp = rbackdir * (-2 * m.dh_b * (rbackmod - m.dh_rc));
			Ftmp.w = m.dh_b * SQR(rbackmod - m.dh_rc);
		}

		// strand ends carry half of the charge
		if(m.dh_half_charged_ends && p_is_end) {
			Ftmp *= (c_number) 0.5;
		}
		if(m.dh_half_charged_ends && q_is_end) {
			Ftmp *= (c_number) 0.5;
		}

		Ttmp -= _cross(ppos_back, Ftmp);
		F -= Ftmp;
	}

	// the angles that enter hydrogen bonding and cross stacking
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	const bool in_hb_range = SQR(DRH_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(DRH_HYDR_RCHIGH);
	const bool in_crst_range = SQR(DRH_CRST_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(DRH_CRST_RCHIGH);
	if(in_hb_range || in_crst_range) {
		c_number rhydromod = sqrt(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		c_number cost1 = -CUDA_DOT(a1, b1);
		c_number cost2 = -CUDA_DOT(b1, rhydrodir);
		c_number cost3 = CUDA_DOT(a1, rhydrodir);
		c_number cost4 = CUDA_DOT(a3, b3);
		c_number cost7 = -CUDA_DOT(b3, rhydrodir);
		c_number cost8 = CUDA_DOT(a3, rhydrodir);

		// HYDROGEN BONDING
		if(pbtype + qbtype == 3 && in_hb_range) {
			const c_number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? m.hb_multiplier : (c_number) 1;
			// as in the CPU code, q plays the role of n3 and p the one of n5 (q is the DNA nucleotide if p is the RNA one)
			const int f1_index = _hybrid_f1_index(qtype, ptype, p_is_rna, !p_is_rna);

			c_number f1 = hb_multi * _hybrid_f1(rhydromod, f1_index);
			c_number f4t1, f4t1D, f4t2, f4t2D, f4t3, f4t3D, f4t4, f4t4D, f4t7, f4t7D, f4t8, f4t8D;
			_f4_of_cos(cost1, HYDR_THETA1, f4t1, f4t1D);
			_f4_of_cos(cost2, HYDR_THETA2, f4t2, f4t2D);
			_f4_of_cos(cost3, HYDR_THETA2, f4t3, f4t3D);
			_f4_of_cos(cost4, HYDR_THETA4, f4t4, f4t4D);
			_f4_of_cos(cost7, HYDR_THETA7, f4t7, f4t7D);
			_f4_of_cos(cost8, HYDR_THETA7, f4t8, f4t8D);

			c_number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
			tot_Tw += hb_energy;

			// the modulations can be exactly 0
			if(hb_energy != (c_number) 0) {
				c_number f1D = hb_multi * _hybrid_f1D(rhydromod, f1_index);
				// the derivatives are the ones of the CPU code (note the signs)
				c_number f4t3Dsin = -f4t3D;
				c_number f4t4Dsin = -f4t4D;
				c_number f4t8Dsin = -f4t8D;

				// RADIAL PART
				c_number4 force = -rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);
				c_number4 torquep = make_c_number4(0, 0, 0, 0);

				// THETA4; t4 = LRACOS (a3 * b3);
				c_number4 dir = _cross(a3, b3);
				c_number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;
				torquep -= dir * torquemod;

				// THETA1; t1 = LRACOS (-a1 * b1);
				dir = _cross(a1, b1);
				torquemod = -f1 * f4t1D * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
				torquep -= dir * torquemod;

				// THETA2; t2 = LRACOS (-b1 * rhydrodir); only q feels a torque
				c_number fact = f1 * f4t1 * f4t2D * f4t3 * f4t4 * f4t7 * f4t8;
				force += (b1 + rhydrodir * cost2) * (fact / rhydromod);

				// THETA3; t3 = LRACOS (a1 * rhydrodir);
				force += (a1 - rhydrodir * cost3) * (f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);
				torquemod = -f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
				torquep += _cross(rhydrodir, a1) * torquemod;

				// THETA7; t7 = LRACOS (-rhydrodir * b3); only q feels a torque
				force += (b3 + rhydrodir * cost7) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7D * f4t8 / rhydromod);

				// THETA8; t8 = LRACOS (rhydrodir * a3);
				force += (a3 - rhydrodir * cost8) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);
				torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
				torquep += _cross(rhydrodir, a3) * torquemod;

				torquep -= _cross(ppos_base, force);

				// the CPU code applies -force to p
				Ftmp = -force;
				Ftmp.w = hb_energy;
				F += Ftmp;
				Ttmp += torquep;
			}
		}

		// CROSS STACKING
		if(in_crst_range) {
			c_number f2 = _hybrid_f2(rhydromod, 0);
			c_number f4t1, f4t1D, f4t2, f4t2D, f4t3, f4t3D, f4t4, f4t4D, f4t7, f4t7D, f4t8, f4t8D;
			_f4_of_cos(cost1, CRST_THETA1, f4t1, f4t1D);
			_f4_of_cos(cost2, CRST_THETA2, f4t2, f4t2D);
			_f4_of_cos(cost3, CRST_THETA2, f4t3, f4t3D);
			_f4_of_cos_symmetrised(cost4, CRST_THETA4, f4t4, f4t4D);
			_f4_of_cos_symmetrised(cost7, CRST_THETA7, f4t7, f4t7D);
			_f4_of_cos_symmetrised(cost8, CRST_THETA7, f4t8, f4t8D);

			c_number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
			if(cstk_energy != (c_number) 0) {
				c_number f2D = _hybrid_f2D(rhydromod, 0);
				c_number f4t3Dsin = -f4t3D;
				c_number f4t4Dsin = -f4t4D;
				c_number f4t8Dsin = -f4t8D;

				// RADIAL PART
				c_number4 force = -rhydrodir * (f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);
				c_number4 torquep = make_c_number4(0, 0, 0, 0);

				// THETA1
				c_number torquemod = -f2 * f4t1D * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
				torquep -= _cross(a1, b1) * torquemod;

				// THETA2: only q feels a torque
				force += (b1 + rhydrodir * cost2) * (f2 * f4t1 * f4t2D * f4t3 * f4t4 * f4t7 * f4t8 / rhydromod);

				// THETA3
				force += (a1 - rhydrodir * cost3) * (f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);
				torquemod = -f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
				torquep += _cross(rhydrodir, a1) * torquemod;

				// THETA4
				torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;
				torquep -= _cross(a3, b3) * torquemod;

				// THETA7: only q feels a torque
				force += (b3 + rhydrodir * cost7) * (f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7D * f4t8 / rhydromod);

				// THETA8
				force += (a3 - rhydrodir * cost8) * (f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);
				torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
				torquep += _cross(rhydrodir, a3) * torquemod;

				torquep -= _cross(ppos_base, force);

				Ftmp = -force;
				Ftmp.w = cstk_energy;
				F += Ftmp;
				Ttmp += torquep;
			}
		}
	}

	// COAXIAL STACKING
	c_number4 rstack = r + qpos_stack - ppos_stack;
	c_number rstackmodsqr = CUDA_DOT(rstack, rstack);
	if(SQR(DRH_CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(DRH_CXST_RCHIGH)) {
		c_number rstackmod = sqrt(rstackmodsqr);
		c_number4 rstackdir = rstack / rstackmod;

		// angles involved in the CXST interaction
		c_number cost5 = CUDA_DOT(a3, rstackdir);
		c_number cost6 = -CUDA_DOT(b3, rstackdir);
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t5 = CUDA_LRACOS(cost5);
		c_number t6 = CUDA_LRACOS(cost6);

		c_number4 rbackbone_cxst = r + qpos_back - ppos_back;
		c_number rbackmod_cxst = _module(rbackbone_cxst);
		c_number4 rbackbonedir = rbackbone_cxst / rbackmod_cxst;
		c_number cosphi3 = CUDA_DOT(rstackdir, _cross(rbackbonedir, a1));
		c_number cosphi4 = CUDA_DOT(rstackdir, _cross(rbackbonedir, b1));

		c_number f2 = _hybrid_f2(rstackmod, 1);
		c_number f4t1 = _f4(t1, CXST_THETA1) + _f4(2 * PI - t1, CXST_THETA1);
		c_number f4t4 = _f4(t4, CXST_THETA4);
		c_number f4t5 = _f4(t5, CXST_THETA5) + _f4(PI - t5, CXST_THETA5);
		c_number f4t6 = _f4(t6, CXST_THETA5) + _f4(PI - t6, CXST_THETA5);
		c_number f5cosphi3 = _f5(cosphi3, DRH_CXST_PHI3_A, DRH_CXST_PHI3_B, DRH_CXST_PHI3_XC, DRH_CXST_PHI3_XS);
		c_number f5cosphi4 = _f5(cosphi4, DRH_CXST_PHI4_A, DRH_CXST_PHI3_B, DRH_CXST_PHI4_XC, DRH_CXST_PHI4_XS);

		c_number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
		if(cxst_energy != (c_number) 0) {
			c_number f2D = _hybrid_f2D(rstackmod, 1);
			c_number f4t1Dsin = -_f4Dsin(t1, CXST_THETA1) + _f4Dsin(2 * PI - t1, CXST_THETA1);
			c_number f4t4Dsin = _f4Dsin(t4, CXST_THETA4);
			c_number f4t5Dsin = _f4Dsin(t5, CXST_THETA5) - _f4Dsin(PI - t5, CXST_THETA5);
			c_number f4t6Dsin = -_f4Dsin(t6, CXST_THETA5) + _f4Dsin(PI - t6, CXST_THETA5);
			c_number f5Dcosphi3 = _f5D(cosphi3, DRH_CXST_PHI3_A, DRH_CXST_PHI3_B, DRH_CXST_PHI3_XC, DRH_CXST_PHI3_XS);
			c_number f5Dcosphi4 = _f5D(cosphi4, DRH_CXST_PHI4_A, DRH_CXST_PHI3_B, DRH_CXST_PHI4_XC, DRH_CXST_PHI4_XS);

			// RADIAL PART
			c_number4 force = -rstackdir * f2D * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
			c_number4 torquep = make_c_number4(0, 0, 0, 0);

			// THETA1; t1 = LRACOS (-a1 * b1);
			c_number torquemod = -f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
			torquep -= _cross(a1, b1) * torquemod;

			// THETA4; t4 = LRACOS (a3 * b3);
			torquemod = -f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
			torquep -= _cross(a3, b3) * torquemod;

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			c_number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * f5cosphi3 * f5cosphi4;
			force += (a3 - rstackdir * cost5) * (fact / rstackmod);
			torquep -= _cross(rstackdir, a3) * fact;

			// THETA6; t6 = LRACOS (-b3 * rstackdir); only q feels a torque
			fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * f5cosphi3 * f5cosphi4;
			force += (b3 + rstackdir * cost6) * (fact / rstackmod);

			torquep -= _cross(ppos_stack, force);

			// cos(phi3) and cos(phi4)
			c_number4 myforce = make_c_number4(0, 0, 0, 0);
			c_number4 mytorque = make_c_number4(0, 0, 0, 0);

			c_number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi4 * f5Dcosphi3;
			c_number4 rba1 = _cross(rbackbonedir, a1);
			c_number4 a1rs = _cross(a1, rstackdir);
			c_number rb_dot_a1rs = CUDA_DOT(rbackbonedir, a1rs);
			c_number rs_dot_rba1 = CUDA_DOT(rstackdir, rba1);
			c_number4 forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
			c_number4 forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod_cxst);
			myforce += forcestack;
			myforce += forceback;
			mytorque -= _cross(ppos_stack, forcestack);
			mytorque -= _cross(ppos_back, forceback);
			mytorque -= _cross(a1, _cross(rstackdir, rbackbonedir)) * force_c;

			force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;
			rba1 = _cross(rbackbonedir, b1);
			a1rs = _cross(b1, rstackdir);
			rb_dot_a1rs = CUDA_DOT(rbackbonedir, a1rs);
			rs_dot_rba1 = CUDA_DOT(rstackdir, rba1);
			forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
			forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod_cxst);
			myforce += forcestack;
			myforce += forceback;
			mytorque -= _cross(ppos_stack, forcestack);
			mytorque -= _cross(ppos_back, forceback);

			force += myforce;
			torquep += mytorque;

			Ftmp = -force;
			Ftmp.w = cxst_energy;
			F += Ftmp;
			Ttmp += torquep;
		}
	}

	T += Ttmp;

	// this component stores the energy due to hydrogen bonding
	T.w = tot_Tw;
}

/**
 * @brief Parameters of the DNA and RNA parts of the interaction that are not stored in constant memory.
 */
struct KernelParameters {
	bool average;
	bool mismatch_repulsion;
	bool dna_use_mbf;
	c_number dna_mbf_xmax, dna_mbf_finf;
	bool rna_use_mbf;
	c_number rna_mbf_xmax, rna_mbf_finf;
	/// testing aid: only the DNA/RNA interactions are computed
	bool hybrid_only;
};

/**
 * @brief Adds to F and T the contribution of the bonded interactions (FENE, bonded excluded volume, stacking) of the nucleotide whose
 * bonds are pbonds. Strands are either DNA or RNA. The energy goes to F.w and the torque is in the frame of the nucleotide.
 */
__device__ void _drh_bonded_part(const c_number4 *poss, const GPU_quat *orientations, const LR_bonds &pbonds, c_number4 &ppos, c_number4 &a1,
		c_number4 &a2, c_number4 &a3, bool p_is_rna, const KernelParameters &params, c_number4 &F, c_number4 &T) {
	if(pbonds.n3 != P_INVALID) {
		c_number4 qpos = poss[pbonds.n3];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n3], b1, b2, b3);

		if(p_is_rna) {
			cuda_rna::_bonded_part<true>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, params.average, params.rna_use_mbf, params.rna_mbf_xmax, params.rna_mbf_finf);
		}
		else {
			c_number4 dF = make_c_number4(0, 0, 0, 0);
			cuda_dna::_bonded_part<true>(qpos - ppos, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, T, true, true, params.dna_use_mbf, params.dna_mbf_xmax, params.dna_mbf_finf);
			F += dF;
		}
	}
	if(pbonds.n5 != P_INVALID) {
		c_number4 qpos = poss[pbonds.n5];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n5], b1, b2, b3);

		if(p_is_rna) {
			cuda_rna::_bonded_part<false>(qpos, b1, b2, b3, ppos, a1, a2, a3, F, T, params.average, params.rna_use_mbf, params.rna_mbf_xmax, params.rna_mbf_finf);
		}
		else {
			c_number4 dF = make_c_number4(0, 0, 0, 0);
			cuda_dna::_bonded_part<false>(ppos - qpos, qpos, b1, b2, b3, ppos, a1, a2, a3, dF, T, true, true, params.dna_use_mbf, params.dna_mbf_xmax, params.dna_mbf_finf);
			F += dF;
		}
	}
}

/**
 * @brief Computes forces and torques for the oxNA interaction using neighbour matrices. The nucleotide with index i is DNA if is_rna[i] == 0 and
 * RNA otherwise, where i is the index that the nucleotide has in the topology (i.e. it does not change if particles are sorted).
 */
__global__ void drh_forces(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs, const int *number_neighs, const LR_bonds __restrict__ *bonds,
		const int * __restrict__ is_rna, KernelParameters params, const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];
	bool p_is_end = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);
	bool p_is_rna = is_rna[get_particle_index(ppos)];

	// particle axes according to Allen's paper
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(!params.hybrid_only) {
		_drh_bonded_part(poss, orientations, pbonds, ppos, a1, a2, a3, p_is_rna, params, F, T);
	}

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);

	T.w = (c_number) 0;
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND && pbonds.n3 != k_index && pbonds.n5 != k_index) {
			c_number4 qpos = poss[k_index];
			c_number4 r = box->minimum_image(ppos, qpos);

			c_number4 b1, b2, b3;
			get_vectors_from_quat(orientations[k_index], b1, b2, b3);
			LR_bonds qbonds = bonds[k_index];
			bool q_is_end = (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID);
			bool q_is_rna = is_rna[get_particle_index(qpos)];

			if(p_is_rna != q_is_rna) {
				_particle_particle_hybrid_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, p_is_rna, p_is_end, q_is_end);
			}
			else if(params.hybrid_only) {
				continue;
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
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

/**
 * @brief Non-bonded part of the oxNA interaction for the edge-based approach: each thread handles one pair of nucleotides of the edge list
 * (which does not contain bonded pairs and contains each pair once). The interaction is evaluated from the point of view of the first nucleotide
 * of the pair, and the force and torque acting on the second one are obtained from the conservation of linear and angular momentum.
 *
 * The torques are accumulated in the lab frame, like the neighbour-matrix kernel does before changing frame: they are converted to the frame of
 * each nucleotide by drh_forces_edge_bonded(), which has to be called afterwards. Unlike the other edge-based kernels, the energies are accumulated too:
 * as in drh_forces, each nucleotide receives in the fourth component of the force the energy of each of its pairs, and in the fourth component of the
 * torque the hydrogen-bonding energy of each of its pairs (each pair is therefore counted twice).
 */
__global__ void drh_forces_edge_nonbonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const edge_bond __restrict__ *edge_list, int n_edges, const LR_bonds __restrict__ *bonds,
		const int * __restrict__ is_rna, KernelParameters params, const CUDABox *box) {
	if(IND >= n_edges) return;

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];
	bool p_is_rna = is_rna[get_particle_index(ppos)];

	// get info for particle 2
	c_number4 qpos = poss[b.to];
	bool q_is_rna = is_rna[get_particle_index(qpos)];

	if(p_is_rna == q_is_rna && params.hybrid_only) return;

	LR_bonds pbonds = bonds[b.from];
	LR_bonds qbonds = bonds[b.to];
	bool p_is_end = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);
	bool q_is_end = (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID);

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[b.from], a1, a2, a3);
	c_number4 b1, b2, b3;
	get_vectors_from_quat(orientations[b.to], b1, b2, b3);

	c_number4 r = box->minimum_image(ppos, qpos);

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 dT = make_c_number4(0, 0, 0, 0);
	if(p_is_rna != q_is_rna) {
		_particle_particle_hybrid_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, p_is_rna, p_is_end, q_is_end);
	}
	else if(p_is_rna) {
		cuda_rna::_particle_particle_RNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, params.average, true, params.mismatch_repulsion, p_is_end, q_is_end);
	}
	else {
		cuda_dna::_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, true, true, true, p_is_end, q_is_end);
	}

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;

	// the fourth components are used to store the energies
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

	// the torque may be different from 0 even if dF == 0
	if(CUDA_DOT(dT, dT) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[to_index]), dT);

	if(pair_energy != (c_number) 0) {
		atomicAdd(&(forces[from_index].w), pair_energy);
		atomicAdd(&(forces[to_index].w), pair_energy);
	}
	if(hb_energy != (c_number) 0) {
		atomicAdd(&(torques[from_index].w), hb_energy);
		atomicAdd(&(torques[to_index].w), hb_energy);
	}
}

/**
 * @brief Bonded part of the oxNA interaction for the edge-based approach. It also changes the frame of the torques from the one of the lab
 * (which is the one in which drh_forces_edge_nonbonded() accumulates them) to the one of each nucleotide.
 */
__global__ void drh_forces_edge_bonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const LR_bonds __restrict__ *bonds, const int * __restrict__ is_rna, KernelParameters params) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];
	bool p_is_rna = is_rna[get_particle_index(ppos)];

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(!params.hybrid_only) {
		_drh_bonded_part(poss, orientations, pbonds, ppos, a1, a2, a3, p_is_rna, params, F, T);
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

} // namespace cuda_drh

} // anonymous namespace

#endif /* CUDA_DRH_CUH_ */

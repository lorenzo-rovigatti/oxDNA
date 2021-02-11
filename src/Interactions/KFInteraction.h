/*
 * KFInteraction.h
 *
 *  Created on: 19/sep/2016
 *      Author: lorenzo
 */

#ifndef KFINTERACTION_H_
#define KFINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Manages the interaction between patchy particles interacting through a KF potential (or a continuous version of it)
 *
 * This interaction is selected with
 * interaction_type = KF
 *
 * @verbatim
 KF_N = <int> (number of patches)
 KF_continuous = <bool> (selects either the original KF model, valid only for MC simulations, or its continuous variant (see ACS Nano 10, 5459 (2016)))
 KF_delta = <float> (radial width of the patches)
 KF_cosmax = <float> (angular half-width of the patches)
 [KF_N_B = <int> (number of patches on species B)]
 [KF_epsilon_AA = <float> (depth of the well of the patch-patch interaction between particles of species A)]
 [KF_epsilon_BB = <float> (depth of the well of the patch-patch interaction between particles of species B)]
 [KF_epsilon_AB = <float> (depth of the well of the patch-patch interaction between particles of unlike species)]
 [KF_sigma_AA = <float> (diameter controlling the repulsive interaction between particles of species A)]
 [KF_sigma_BB = <float> (diameter controlling the repulsive interaction between particles of species B)]
 [KF_sigma_AB = <float> (diameter controlling the repulsive interaction between particles of unlike species)]
 @endverbatim
 */

class KFInteraction: public BaseInteraction {
protected:
	/// Number of patches per particle
	int _N_patches;

	/// Number of patches per second-species particle
	int _N_patches_B;

	/// Number of particles per species
	int _N_A, _N_B;

	/// True if we are to simulate a patchy binary mixture, false otherwise
	bool _is_binary;

	/// Trus if we are to simulate non-hard bodies
	bool _is_continuous;

	/// Exponent of the inverse-power-law repulsion between the particles
	int _rep_power;

	/// Particles' diameters
	number _sigma[3];

	/// Squared diameters
	number _sqr_sigma[3];

	/// Repulsive interaction cut-off squared
	number _sqr_tot_rcut[3];

	/// Depth of the patch-patch well
	number _epsilon[3];

	/// Repulsive interaction energy at the cut-off
	number _E_cut[3];

	/// Exponent for the Gaussian-like potential well used for the patches
	int _patch_power;

	/// Width of the patches
	number _patch_delta;

	/// Angular width of the patches
	number _patch_cosmax;

	/// _patch_alpha^10
	number _patch_pow_delta;

	/// _patch_cosmax^30
	number _patch_pow_cosmax;

	/// Angular cut-off for the patchy attraction
	number _patch_angular_cutoff;

	/**
	 * @brief KF interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	inline number _KF_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	/**
	 * @brief Continuous KF interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	inline number _continuous_KF_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

public:
	enum {
		PATCHY = 4
	};

	KFInteraction();
	virtual ~KFInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

number KFInteraction::_continuous_KF_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	int type = p->type + q->type;
	if(sqr_r > _sqr_tot_rcut[type])
		return (number) 0.f;

	number energy = (number) 0.f;

	number part = pow(_sqr_sigma[type] / sqr_r, _rep_power * 0.5);
	energy = part - _E_cut[type];

	if(update_forces) {
		LR_vector force = _computed_r * (_rep_power * part / sqr_r);
		p->force -= force;
		q->force += force;
	}

	// here everything is done as in Allen's paper
	number rmod = sqrt(sqr_r);
	LR_vector r_versor = _computed_r / (-rmod);

	number sqr_surf_dist = SQR(rmod - 1.);
	number r8b10 = SQR(SQR(sqr_surf_dist)) / _patch_pow_delta;
	number exp_part = -1.001 * exp(-(number) 0.5 * r8b10 * sqr_surf_dist);

	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
		LR_vector ppatch = p->int_centers[pi] * 2.;

		number cospr = -(ppatch * r_versor);
		if(cospr > _patch_angular_cutoff) {
			number cospr_base = pow(cospr - 1., _patch_power - 1);
			// we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
			number cospr_part = cospr_base * (cospr - 1.);
			number p_mod = exp(-cospr_part / (2. * _patch_pow_cosmax));

			for(uint pj = 0; pj < q->N_int_centers(); pj++) {
				LR_vector qpatch = q->int_centers[pj] * 2.;

				number cosqr = qpatch * r_versor;
				if(cosqr > _patch_angular_cutoff) {
					number cosqr_base = pow(cosqr - 1., _patch_power - 1);
					number cosqr_part = cosqr_base * (cosqr - 1.);
					number q_mod = exp(-cosqr_part / (2. * _patch_pow_cosmax));

					energy += exp_part * p_mod * q_mod;

					if(update_forces) {
						// radial part
						LR_vector tmp_force = r_versor * (p_mod * q_mod * 5. * (rmod - 1.) * exp_part * r8b10);

						// angular p part
						number der_p = exp_part * q_mod * (0.5 * _patch_power * p_mod * cospr_base / _patch_pow_cosmax);
						LR_vector p_ortho = ppatch + cospr * r_versor;
						tmp_force -= p_ortho * (der_p / rmod);

						// angular q part
						number der_q = exp_part * p_mod * (-0.5 * _patch_power * q_mod * cosqr_base / _patch_pow_cosmax);
						LR_vector q_ortho = qpatch - cosqr * r_versor;
						tmp_force -= q_ortho * (der_q / rmod);

						p->force += tmp_force;
						q->force -= tmp_force;

						p->torque += p->orientationT * (r_versor.cross(ppatch) * der_p);
						q->torque += q->orientationT * (r_versor.cross(qpatch) * der_q);
					}
				}
			}
		}
	}

	return energy;
}

number KFInteraction::_KF_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(update_forces)
		throw oxDNAException("KFInteraction: forces are not defined in non-continuous KF interactions");

	number sqr_r = _computed_r.norm();
	int type = p->type + q->type;

	if(sqr_r > _sqr_tot_rcut[type]) {
		return 0.f;
	}
	if(sqr_r < _sqr_sigma[type]) {
		this->set_is_infinite(true);
		return 1e10;
	}

	number energy = (number) 0.f;
	LR_vector r_versor(_computed_r / sqrt(sqr_r));
	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
		// the factor of two comes from the fact that patches are normalised to 0.5
		LR_vector ppatch = p->int_centers[pi] * 2.;

		number p_cos = r_versor * ppatch;
		if(p_cos > _patch_cosmax) {
			for(uint pj = 0; pj < q->N_int_centers(); pj++) {
				LR_vector qpatch = q->int_centers[pj] * 2.;

				number q_cos = -r_versor * qpatch;
				if(q_cos > _patch_cosmax) {
					energy -= _epsilon[type];
				}
			}
		}
	}

	return energy;
}

#endif /* KFINTERACTION_H_ */

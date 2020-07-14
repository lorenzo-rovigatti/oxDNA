/*
 * NathanStarInteraction.h
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#ifndef NATHANSTARINTERACTION_H_
#define NATHANSTARINTERACTION_H_

#include "Interactions/BaseInteraction.h"

#include <gsl/gsl_spline.h>

class NathanStarInteraction: public BaseInteraction<NathanStarInteraction> {
protected:
	number _T;
	const number _sqrt_pi;

	/// Repulsive interaction energy at the cut-off
	number _rep_E_cut;
	int _rep_power;

	/// Exponent for the Guassian-like potential well used for the patches
	int _patch_power;
	/// Cosine of the patch half width
	number _patch_cosmax;
	/// Cosmax to the power of _patch_power
	number _patch_pow_sigma;
	/// Cut-off for the attractive part. It is computed at runtime
	number _patch_cutoff;
	/// Angular cut-off for the patchy attraction
	number _patch_angular_cutoff;
	/// Width of the patch, defaults to 0.12
	number _patch_alpha;
	/// _patch_alpha^10
	number _patch_pow_alpha;
	/// particle radius
	const number _patch_r, _sqr_patch_r;
	/// Interaction cut-off between patchy particles
	number _patchy_rcut, _sqr_patchy_rcut;

	number _star_sigma_g, _sqr_star_sigma_g, _star_sigma_s, _sqr_star_sigma_s;
	number _star_rg, _star_rs, _sqr_star_rs;
	number _star_factor;
	/// Interaction cut-off between patchy particles and stars
	number _patchy_star_rcut, _sqr_patchy_star_rcut;
	/// Interaction cut-off between stars
	number _star_rcut, _sqr_star_rcut;
	int _star_f;
	number _star_f1_2, _star_f3_2;

	number _pi_lambda;
	number _kappa_rs, _kappa;
	number _xi, _zeta;
	number _exp_sqr_kappa_rs;
	number _erfc_kappa_rs;

	int _N_patchy;
	int _N_stars;
	bool _is_marzi;
	bool _make_crystal;
	int _N_in_crystal;
	std::string _crystal_type;

	void _setup_lambda_kappa();

	number _psi_1(number x, number smax);
	number _psi_2(number x, number smax);

	int _interp_size;
	gsl_spline *_spl_patchy, *_spl_patchy_star;
	gsl_interp_accel *_acc_patchy, *_acc_patchy_star;

	number _pressure(number s);
	number _patchy_star_derivative(number r);
	number _patchy_star_marzi_derivative(number, gsl_spline *, gsl_interp_accel *, number);
	void _setup_interp();

	number _patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _patchy_star_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _star_star_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
public:
	enum {
		PATCHY_PATCHY = 0, PATCHY_POLYMER = 1, POLYMER_POLYMER = 2
	};

	enum {
		PATCHY_PARTICLE = 0, POLYMER = 1
	};

	NathanStarInteraction();
	virtual ~NathanStarInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);
	virtual int get_N_from_topology();
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles);
};

class NathanPatchyParticle: public BaseParticle {
public:
	NathanPatchyParticle() :
					BaseParticle() {
	}
	;
	virtual ~NathanPatchyParticle() {
	}
	;

	virtual bool is_rigid_body() {
		return true;
	}
};

class NathanPolymerParticle: public BaseParticle {
public:
	NathanPolymerParticle() :
					BaseParticle() {
	}
	;
	virtual ~NathanPolymerParticle() {
	}
	;

	virtual bool is_rigid_body() {
		return false;
	}
	virtual bool is_bonded(BaseParticle *q) {
		if(q->type == NathanStarInteraction::POLYMER) {
			return (q == this->n3 || q == this->n5);
		}
		return false;
	}
};

extern "C" NathanStarInteraction *make_NathanStarInteraction();

#endif /* NATHANSTARINTERACTION_H_ */

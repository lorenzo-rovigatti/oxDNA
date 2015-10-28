/*
 * NathanStarInteraction.h
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#ifndef NATHANSTARINTERACTION_H_
#define NATHANSTARINTERACTION_H_

#include "BaseInteraction.h"

#include <gsl/gsl_spline.h>

template <typename number>
class NathanStarInteraction: public BaseInteraction<number, NathanStarInteraction<number> > {
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
	string _crystal_type;

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

	number _patchy_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _patchy_star_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _star_star_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
public:
	enum {
		PATCHY_PATCHY = 0,
		PATCHY_POLYMER = 1,
		POLYMER_POLYMER = 2
	};

	enum {
		PATCHY_PARTICLE = 0,
		POLYMER = 1
	};

	NathanStarInteraction();
	virtual ~NathanStarInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual int get_N_from_topology();
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
	virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

template<typename number>
number NathanStarInteraction<number>::_patchy_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->type != PATCHY_PARTICLE && q->type != PATCHY_PARTICLE) return 0.f;
	number rnorm = r->norm();
	if(rnorm > _sqr_patchy_rcut) return (number) 0.f;

	// here everything is done as in Allen's paper
	number rmod = sqrt(rnorm);
	LR_vector<number> r_versor = *r / (-rmod);

	// repulsion
	number rep_part = 1. / pow(rnorm, _rep_power / 2);
	number energy = rep_part - _rep_E_cut;
	if(update_forces) {
		LR_vector<number> force = r_versor * (_rep_power * rep_part / rmod);
		p->force += force;
		q->force -= force;
	}

	// attraction
	LR_vector<number> p_axis = p->orientationT.v3;
	LR_vector<number> q_axis = q->orientationT.v3;

	number cospr = -(p_axis*r_versor);
	if(cospr < 0.) {
		p_axis = -p_axis;
		cospr = -cospr;
	}
	number cosqr = q_axis*r_versor;
	if(cosqr < 0.) {
		q_axis = -q_axis;
		cosqr = -cosqr;
	}
	if(cospr < _patch_angular_cutoff || cosqr < _patch_angular_cutoff) return energy;

	number cospr_base = pow(cospr - 1., _patch_power - 1);
	// we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
	number cospr_part = cospr_base*(cospr - 1.);
	number p_mod = exp(-cospr_part / (2.*_patch_pow_sigma));

	number cosqr_base = pow(cosqr - 1., _patch_power - 1);
	number cosqr_part = cosqr_base*(cosqr - 1.);
	number q_mod = exp(-cosqr_part / (2.*_patch_pow_sigma));

	number sqr_surf_dist = SQR(rmod - 1.);
	number r8b10 = SQR(SQR(sqr_surf_dist)) / _patch_pow_alpha;
	number exp_part = -1.001 * exp(-(number)0.5 * r8b10 * sqr_surf_dist);
	energy += exp_part*p_mod*q_mod;

	if(update_forces) {
		// radial part
		LR_vector<number> tmp_force = r_versor * (p_mod*q_mod * 5.*(rmod - 1.)*exp_part*r8b10);

		// angular p part
		number der_p = exp_part * q_mod * (0.5*_patch_power*p_mod * cospr_base / _patch_pow_sigma);
		LR_vector<number> p_ortho = p_axis + cospr*r_versor;
		tmp_force -= p_ortho * (der_p/rmod);

		// angular q part
		cosqr_part /= cosqr - 1.;
		number der_q = exp_part * p_mod * (-0.5*_patch_power*q_mod * cosqr_base / _patch_pow_sigma);
		LR_vector<number> q_ortho = q_axis - cosqr*r_versor;
		tmp_force -= q_ortho * (der_q/rmod);

		p->force += tmp_force;
		q->force -= tmp_force;

		p->torque += p->orientationT * (r_versor.cross(p_axis) * der_p);
		q->torque += q->orientationT * (r_versor.cross(q_axis) * der_q);
	}

	return energy;
}

template<typename number>
number NathanStarInteraction<number>::_patchy_star_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = (number) 0.f;
	number sqr_r = r->norm();
	if(sqr_r > _sqr_patchy_star_rcut) return (number) 0.f;

	number mod_r = sqrt(sqr_r);

	// this is just to avoid spitting out NaNs and it should occur very rarely, and only during equilibration
	if(mod_r < _spl_patchy_star->interp->xmin) energy = 1e6;
	// this can happen only in single precision
	else if(mod_r > _spl_patchy_star->interp->xmax) return (number) 0.f;
	else energy = gsl_spline_eval(_spl_patchy_star, mod_r, _acc_patchy_star);
	if(update_forces) {
		number force_mod = (mod_r < _spl_patchy_star->interp->xmin) ? 100 : -gsl_spline_eval_deriv(_spl_patchy_star, mod_r, _acc_patchy_star)/mod_r;
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}
	return energy;
}

template<typename number>
number NathanStarInteraction<number>::_star_star_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = (number) 0.f;
	number sqr_r = r->norm();
	if(sqr_r > _sqr_star_rcut) return (number) 0.f;

	number mod_r = sqrt(sqr_r);

	number common_fact = _star_factor * 5. * _T * _star_f3_2 / 18.;
	number i_f = 1. / (1. + _star_f1_2*0.5);

	if(sqr_r < _sqr_star_sigma_s) {
		energy = -log(mod_r / _star_sigma_s) + i_f;

		if(update_forces) {
			// force over r
			number force_mod = common_fact / sqr_r;
			p->force -= *r * force_mod;
			q->force += *r * force_mod;
		}
	}
	else {
		number exp_factor = exp(-(mod_r - _star_sigma_s)*_star_f1_2/(2.*_star_sigma_s));
		energy = i_f * (_star_sigma_s / mod_r) * exp_factor;

		if(update_forces) {
			// force over r
			number force_mod = common_fact * i_f * exp_factor * (_star_sigma_s/(sqr_r*mod_r) + _star_f1_2/(2.*sqr_r));
			p->force -= *r * force_mod;
			q->force += *r * force_mod;
		}
	}

	return common_fact*energy;
}

template<typename number>
class NathanPatchyParticle: public BaseParticle<number> {
public:
	NathanPatchyParticle() : BaseParticle<number>() {};
	virtual ~NathanPatchyParticle() {};

	virtual bool is_rigid_body() { return true; }
};

template<typename number>
class NathanPolymerParticle: public BaseParticle<number> {
public:
	NathanPolymerParticle() : BaseParticle<number>() {};
	virtual ~NathanPolymerParticle() {};

	virtual bool is_rigid_body() { return false; }
	virtual bool is_bonded(BaseParticle<number> *q) {
		if(q->type == NathanStarInteraction<number>::POLYMER) {
			return (q == this->n3 || q == this->n5);
		}
		return false;
	}
};

extern "C" NathanStarInteraction<float> *make_interaction_float();
extern "C" NathanStarInteraction<double> *make_interaction_double();

#endif /* NATHANSTARINTERACTION_H_ */

// OLD METHODS

//template<typename number>
//number NathanStarInteraction<number>::_patchy_star_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//	number energy = (number) 0.f;
//	number sqr_r = r->norm();
//	if(sqr_r > _sqr_patchy_star_rcut) return (number) 0.f;
//
//	number mod_r = sqrt(sqr_r);
//	number z = mod_r - _patch_r;
//
//	number common_fact = _T * _pi_lambda * _star_f3_2 * _patch_r / mod_r;
//
//	if(z < _star_rs) {
//		energy = common_fact * (-log(z/_star_rs) - (SQR(z)/_sqr_star_rs - 1.)*(_xi - 0.5) + _zeta);
//		if(z <= 0.) energy = 1e6;
//
//		if(update_forces) {
//			// force over r
//			number smax = sqrt(z*(z + 1.));
//			number force_mod = common_fact / sqr_r;
//			force_mod *= (sqr_r - _sqr_patch_r)*(0.5/SQR(z) - 0.5/_sqr_star_rs + _psi_1(_star_rs, smax)) - log(z/_star_rs) + _psi_2(_star_rs, smax);
//			if(z <= 0.) force_mod = 1000;
//
//			p->force -= *r * force_mod;
//			q->force += *r * force_mod;
//		}
//	}
//	else {
//		energy = common_fact * _zeta * erfc(_kappa*z) / _erfc_kappa_rs;
//
//		if(update_forces) {
//			// force over r
//			number smax = sqrt(z*(z + 1.));
//			number force_mod = common_fact / sqr_r;
//			force_mod *= (sqr_r - _sqr_patch_r)*_psi_1(z, smax) + _psi_2(z, smax);
//
//			p->force -= *r * force_mod;
//			q->force += *r * force_mod;
//		}
//	}
//
//	return energy;
//}

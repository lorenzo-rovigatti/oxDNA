#ifndef TEP_INTERACTION_H
#define TEP_INTERACTION_H

#include "BaseInteraction.h"
#include "../Particles/TEPParticle.h"

/**
 * @brief Handles interactions between elements of the Twistable Elastic Polymer (TEP) model.
 *
 * The model is described in the paper (add reference), eq.4.
 * This description should be edited when the model has finally been developed.
 *
 * Input options:
 *
 * @verbatim
 [_prefer_harmonic_over_fene = <bool>(if True, neighbouring beads are bound by an harmonic potential instead of a FENE one. Defaults to false.)]
 [_allow_broken_fene = <bool>(if True, the code won't die when the beads are beyond the acceptable FENE range. Defaults to True.)]
 @endverbatim
 The option above are true for the DNAinteraction class, not this one, but I keep them so that
 I can keep track of how they should be declared.
 */
template<typename number>
class TEPInteraction: public BaseInteraction<number, TEPInteraction<number> > {
protected:
	number _T;
	//parameters of the TEP model
	number _ka;
	number _kb;
	number _kt;

	number * _kt_pref;
	number * _kb1_pref;
	// parameters of the modified twisting terms
	number _twist_a;
	number _twist_b;

	// parameters of the double bending term

	number * _xu_bending;
	number * _xk_bending;
	number * _kb2_pref;
	number _xu_bending_default;
	number _xk_bending_default;
	number _kb2_pref_default;
	// FENE parameters
	number _TEP_FENE_DELTA;
	number _TEP_FENE_DELTA2;
	number _TEP_FENE_EPS;
	number _TEP_FENE_R0;
	// LJ parameters
	number _TEP_EXCL_EPS_BONDED;
	number _TEP_EXCL_EPS_NONBONDED;
	number _TEP_EXCL_S2;
	number _TEP_EXCL_R2;
	number _TEP_EXCL_B2;
	number _TEP_EXCL_RC2;

	number _TEP_spring_offset;
	/// true by default; set this to false if you want the code to not die when bonded backbones are found to be outside the acceptable FENE range
	bool _allow_broken_fene;

// false by default: set this to true if you want neighbouring particles to be bound by a quadratic potential instead of a FENE
	bool _prefer_harmonic_over_fene;
// zero by default: if nonzero, the beads with a virtual neighbour will be twisted (e.g. will try and align the potential v2 vector to a spinning vector). This twist means that the twist_extremal_particles interaction will be used.
	number _twist_boundary_stiff;

	LR_vector<number> _o1;
	LR_vector<number> _o2;
	number _o1_modulus;
	number _o2_modulus;
	LR_vector<number> _w1;
	LR_vector<number> _w2;
	llint _my_time1;
	llint _my_time2;
	llint _max_twisting_time;
	llint _print_torques_every;

	input_file * input_file_copy;

	virtual number _spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_double_bending(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_twist(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_alignment(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _nonbonded_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	inline number _repulsive_lj2(number prefactor, const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces);
	int setNonNegativeNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description);
	int setNonNegativeLLInt(input_file *inp, const char * skey, llint *dest, int mandatory, const char * arg_description);
	int setPositiveNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description);
	int setNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description);
	int setPositiveLLInt(input_file *inp, const char * skey, llint *dest, int mandatory, const char * arg_description);
	// Return the parameters for the smoothed potential
	// parameter A
	inline number _get_A_smooth(int particle_index) {
		number kb2 = _kb2_pref[particle_index];
		number gu = cos(_xu_bending[particle_index]);
		number gu2 = gu * gu;
		number gu3 = gu2 * gu;
		number gu4 = gu2 * gu2;
		number gk = cos(_xk_bending[particle_index]);
		number gk2 = gk * gk;
		number gk3 = gk2 * gk;
		number gk4 = gk2 * gk2;

		return (gk2 * kb2 - 2. * gk * gu * kb2 + 1. * gu2 * kb2) / (gk4 - 4. * gk3 * gu + 6. * gk2 * gu2 - 4. * gk * gu3 + gu4);
	}

	// parameter B
	inline number _get_B_smooth(int particle_index) {
		number kb2 = _kb2_pref[particle_index];
		number kb1 = _kb1_pref[particle_index];
		number gu = cos(_xu_bending[particle_index]);
		number gu2 = gu * gu;
		number gu3 = gu2 * gu;
		number gu4 = gu2 * gu2;
		number gu5 = gu4 * gu;
		number gk = cos(_xk_bending[particle_index]);
		number gk2 = gk * gk;
		number gk3 = gk2 * gk;
		number gk4 = gk2 * gk2;
		number gk5 = gk4 * gk;
		return (gk * gu3 * (2. * kb1 - 5. * kb2) + gk4 * (-0.5 * kb1 - kb2) + gk3 * gu * (2. * kb1 + kb2) + gu4 * (-0.5 * kb1 + 2. * kb2) + gk2 * gu2 * (-3. * kb1 + 3. * kb2)) / (gk5 - 5. * gk4 * gu + 10. * gk3 * gu2 - 10. * gk2 * gu3 + 5. * gk * gu4 - gu5);

		return (gk * SQR(gu) * gu * (2. * kb1 - 5. * kb2) + SQR(gk) * SQR(gk) * (-0.5 * kb1 - kb2) + SQR(gk) * gk * gu * (2. * kb1 + kb2) + SQR(gu) * SQR(gu) * (-0.5 * kb1 + 2. * kb2) + SQR(gk) * SQR(gu) * (-3. * kb1 + 3. * kb2)) / (SQR(gk) * SQR(gk) * gk - 5. * SQR(gk) * SQR(gk) * gu + 10. * SQR(gk) * gk * SQR(gu) - 10. * SQR(gk) * SQR(gu) * gu + 5. * gk * SQR(gu) * SQR(gu) - SQR(gu) * SQR(gu) * gu);
	}
	// parameter C
	inline number _get_C_smooth(int particle_index) {
		number gu = cos(_xu_bending[particle_index]);
		number gk = cos(_xk_bending[particle_index]);
		number kb2 = _kb2_pref[particle_index];
		number kb1 = _kb1_pref[particle_index];

		return (SQR(gk) * SQR(gk) * gk * kb1 + SQR(gk) * gk * SQR(gu) * (6. * kb1 - 5. * kb2) - SQR(gu) * SQR(gu) * gu * kb2 + gk * SQR(gu) * SQR(gu) * (kb1 + kb2) + SQR(gk) * SQR(gk) * gu * (-4. * kb1 + 2. * kb2) + SQR(gk) * SQR(gu) * gu * (-4. * kb1 + 3. * kb2)) / (SQR(gk) * SQR(gk) * gk - 5. * SQR(gk) * SQR(gk) * gu + 10. * SQR(gk) * gk * SQR(gu) - 10. * SQR(gk) * SQR(gu) * gu + 5. * gk * SQR(gu) * SQR(gu) - SQR(gu) * SQR(gu) * gu);
	}

	// parameter D
	inline number _get_D_smooth(int particle_index) {
		number gu = cos(_xu_bending[particle_index]);
		number kb1 = _kb1_pref[particle_index];

		return kb1 * (1. - gu);
	}
	// parameter phi
	inline number _get_phi_bending(int particle_index) {
		number kb1 = _kb1_pref[particle_index];
		number kb2 = _kb2_pref[particle_index];
		number gu = cos(_xu_bending[particle_index]);
		number gk = cos(_xk_bending[particle_index]);

		return kb1 * (1 - (gu + gk) * 0.5) - kb2 * (1 - gk);
	}
	/**
	 * @brief Checks whether the two particles share a backbone link.
	 *
	 * @param p
	 * @param q
	 * @return true if they are bonded, false otherwise
	 */
	inline bool _are_bonded(BaseParticle<number> *p, BaseParticle<number> *q) {
		return (p->n3 == q || p->n5 == q);
	}

	// The following are all the functions used to implement the twisting of the extremal beads. Hopefully to be removed when torques are programmed properly.
	int getInputDirection(input_file *inp, const char * skey, LR_vector<number> *dest, int mandatory);
	LR_vector<number> rotateVectorAroundVersor(const LR_vector<number> vector, const LR_vector<number> versor, const number angle);
	number _twist_boundary_particles(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _index_twist_boundary_particles(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
public:
	enum {
		SPRING = 0, BONDED_BENDING = 1, BONDED_TWIST = 2, BONDED_ALIGNMENT = 3, NONBONDED_EXCLUDED_VOLUME = 4, BONDED_DEBYE_HUCKEL = 5, NONBONDED_DEBYE_HUCKEL = 6
	};
	TEPInteraction();
	virtual ~TEPInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r = NULL, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r = NULL, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r = NULL, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r = NULL, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	virtual void read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles);

};

template<typename number>
number TEPInteraction<number>::_repulsive_lj2(number prefactor, const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = prefactor * b * SQR(rrc);
			if(update_forces)
				force = -r * (2 * prefactor * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			// the additive term was added by me to mimick Davide's implementation
			energy = 4 * prefactor * (SQR(lj_part) - lj_part) + prefactor;
			if(update_forces)
				force = -r * (24 * prefactor * (lj_part - 2 * SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0)
		force.x = force.y = force.z = (number) 0;

	return energy;
}

template<typename number>
int TEPInteraction<number>::setNonNegativeNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description) {

	int ret_value = getInputNumber(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest < 0)
			throw oxDNAException("read negative parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
	return ret_value;
}

template<typename number>
int TEPInteraction<number>::setPositiveNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description) {
	int ret_value = getInputNumber(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest <= 0)
			throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
	return ret_value;
}

template<typename number>
int TEPInteraction<number>::setNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description) {
	int ret_value = getInputNumber(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
	return ret_value;
}

template<typename number>
int TEPInteraction<number>::setPositiveLLInt(input_file *inp, const char * skey, llint *dest, int mandatory, const char * arg_description) {
	int ret_value = getInputLLInt(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest <= 0)
			throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %lld",skey,*dest);
	}
	return ret_value;
}
template<typename number>
int TEPInteraction<number>::setNonNegativeLLInt(input_file *inp, const char * skey, llint *dest, int mandatory, const char * arg_description) {
	int ret_value = getInputLLInt(inp, skey, dest, mandatory) == KEY_FOUND;
	if(ret_value) {
		if(*dest < 0)
			throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting", skey, arg_description, skey, *dest);
		OX_LOG(Logger::LOG_INFO,"%s manually set to %lld",skey,*dest);
	}
	return ret_value;
}
// The following are all the functions used to implement the twisting of the extremal beads. Hopefully to be removed when torques are programmed properly.
//
// Read a direction from file (TODO: propose to add this to defs.h, since it might be needed elsewhere).
template<typename number>
int TEPInteraction<number>::getInputDirection(input_file *inp, const char * skey, LR_vector<number> *dest, int mandatory) {
	std::string strdir;
	int ret_value = getInputString(inp, skey, strdir, mandatory);
	double x, y, z;
	int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", &x, &y, &z);
	if(tmpi != 3)
		throw oxDNAException("could not parse direction %s in input file. Dying badly", skey);
	*dest = LR_vector<number>((number) x, (number) y, number(z));
	if(x == 0 && y == 0 && z == 0) {
		throw oxDNAException("direction %s in input file is the zero vector.", skey);
	}
	dest->normalize();
	return ret_value;
}

// Rotate a vector around a versor (TODO: propose to add this to defs.h)
template<typename number>
LR_vector<number> TEPInteraction<number>::rotateVectorAroundVersor(const LR_vector<number> vector, const LR_vector<number> versor, const number angle) {
	/* According to section 5.2 of this webpage http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ ,
	 // the image of rotating a vector (x,y,z) around a versor (u,v,w) is given by
	 //      / u(ux + vy + wz)(1 - cos th) + x cos th + (- wy + vz) sin th \
//      | v(ux + vy + wz)(1 - cos th) + y cos th + (+ wx - uz) sin th |
	 //      \ w(ux + vy + wz)(1 - cos th) + z cos th + (- vx + uy) sin th /
	 //      Note that the first parenthesis in every component contains the scalar product of the vectors,
	 //      and the last parenthesys in every component contains a component of the cross product versor ^ vector.
	 //      The derivation that they show on the webpage seems sound, and I've checked it with Mathematica in about//			 1 hour of ininterrupting swearing. */
	number costh = cos(angle);
	number sinth = sin(angle);
	number scalar = vector * versor;
	LR_vector<number> cross = versor.cross(vector);
	return LR_vector<number>(versor.x * scalar * (1. - costh) + vector.x * costh + cross.x * sinth, versor.y * scalar * (1. - costh) + vector.y * costh + cross.y * sinth, versor.z * scalar * (1. - costh) + vector.z * costh + cross.z * sinth);
}
#endif /* TEP_INTERACTION_H */

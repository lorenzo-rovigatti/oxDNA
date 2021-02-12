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

class TEPInteraction: public BaseInteraction {
protected:
	bool _is_on_cuda;

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

	number * _th_b_0;
	number * _beta_0;
	number * _xu_bending;
	number * _xk_bending;
	number * _kb2_pref;

	number _th_b_0_default;
	number _beta_0_default;
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

	LR_vector _o1;
	LR_vector _o2;
	number _o1_modulus;
	number _o2_modulus;
	LR_vector _w1;
	LR_vector _w2;
	int _time_increment;
	llint _my_time1;
	llint _my_time2;
	llint _time_var;
	llint _max_twisting_time;
	llint _choose_direction_time;
	llint _print_torques_every;

	virtual number _spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_bending(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_double_bending(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_twist(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_alignment(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _nonbonded_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	inline number _repulsive_lj2(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);
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
		//return (gk*SQR(gu)*gu*(2.*kb1 - 5.*kb2) + SQR(gk)*SQR(gk)*(-0.5*kb1 - kb2) + SQR(gk)*gk*gu*(2.*kb1 + kb2) + SQR(gu)*SQR(gu)*(-0.5*kb1 + 2.*kb2) +  SQR(gk)*SQR(gu)*(-3.*kb1 + 3.*kb2)) / (SQR(gk)*SQR(gk)*gk - 5.*SQR(gk)*SQR(gk)*gu + 10.*SQR(gk)*gk*SQR(gu) - 10.*SQR(gk)*SQR(gu)*gu + 5.*gk*SQR(gu)*SQR(gu) - SQR(gu)*SQR(gu)*gu);
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
	inline bool _are_bonded(BaseParticle *p, BaseParticle *q) {
		return (p->n3 == q || p->n5 == q);
	}
	/**
	 * @brief update the time_variable used in index_twist_boundary_particles
	 */
	inline llint update_time_variable(llint time) {
		// if we're twisting in a specified direction for a specified time, update it until you reach it
		if(_max_twisting_time != 0) {
			return std::min(time + 1, _max_twisting_time);
		}
		// otherwise, increase or decrease it as needed
		else if(_choose_direction_time != 0) {
			return time + _time_increment;
		}
		// otherwise something don't update it
		return time;
	}

	inline void update_increment(llint time) {
		if(_choose_direction_time != 0) {
			if(time % _choose_direction_time == 0) {
				// time increment can be -1 or 1 with probability 0.5
				_time_increment = drand48() > 0.5 ? -1 : 1;
				FILE * fp = fopen("twist_diff.txt", "a");
				fprintf(fp, "%lld\t%d\n", time, _time_increment);
				fclose(fp);
				//OX_LOG(Logger::LOG_INFO," Setting time increment to %d",_time_increment);
			}
		}
	}

	// The following are all the functions used to implement the twisting of the extremal beads. Hopefully to be removed when torques are programmed properly.
	int getInputDirection(input_file *inp, const char * skey, LR_vector *dest, int mandatory);
	LR_vector rotateVectorAroundVersor(const LR_vector vector, const LR_vector versor, const number angle);
	number _twist_boundary_particles(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _index_twist_boundary_particles(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
public:
	enum {
		SPRING = 0, BONDED_BENDING = 1, BONDED_TWIST = 2, BONDED_ALIGNMENT = 3, NONBONDED_EXCLUDED_VOLUME = 4, BONDED_DEBYE_HUCKEL = 5, NONBONDED_DEBYE_HUCKEL = 6
	};
	TEPInteraction();
	virtual ~TEPInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

};

#endif /* TEP_INTERACTION_H */

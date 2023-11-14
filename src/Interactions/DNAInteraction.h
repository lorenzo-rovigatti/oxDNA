
#ifndef DNA_INTERACTION_H
#define DNA_INTERACTION_H

#include "BaseInteraction.h"

/**
 * @brief Handles interactions between DNA nucleotides.
 *
 * It supports both the average-sequence model (http://jcp.aip.org/resource/1/jcpsa6/v134/i8/p085101_s1)
 * and the sequence-dependence model (http://jcp.aip.org/resource/1/jcpsa6/v137/i13/p135101_s1).
 * This is the default interaction.
 */

class DNAInteraction : virtual public BaseInteraction {
protected:
	bool _average;
	std::string _seq_filename;
	number _T;
	number _hb_multiplier;
	bool _grooving;
	/// true by default; set this to false if you want the code to not die when bonded backbones are found to be outside the acceptable FENE range
	bool _allow_broken_fene;

	/// here we store the r0, which is different in oxDNA and oxDNA2, and is set in the constructor
	number _fene_r0;

	/// variables used when max_backbone_force = true
	bool _use_mbf;
	number _mbf_xmax;
	number _mbf_fmax;
	number _mbf_finf;

	int MESH_F4_POINTS[13];
	Mesh _mesh_f4[13];

	virtual void _on_T_update();

	number _f1(number r, int type, int n3, int n5);
	number _f1D(number r, int type, int n3, int n5);
	number _f2(number r, int type);
	number _f2D(number r, int type);
	number _f4(number t, int type);
	number _f4D(number t, int type);
	number _f4Dsin(number t, int type);
	number _f5(number f, int type);
	number _f5D(number f, int type);

	number _fakef4(number t, void * par);
	number _fakef4D(number t, void * par);
	number _fakef4_cxst_t1(number t, void * par);
	number _fakef4D_cxst_t1(number t, void * par);

	virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	/**
	 * @brief Custom function that returns f4. This was added to add the possibility to avoid the use of meshes in classes that inherit from this.
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4 (number cost, int i) { return this->_mesh_f4[i].query(cost); }

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4D (number cost, int i) { return this->_mesh_f4[i].query_derivative(cost); }

	inline number _repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);

	/**
	 * @brief Check the relation between p and q. Used by the bonded interaction terms.
	 *
	 * Check whether q is the 3' neighbour of p. If this is not the case, it changes identities
	 * between q and p if p is the 3' neighbour of q and updates r accordingly. It return false
	 * if the two particles are not bonded neighbours, true otherwise.
	 * @param p
	 * @param q
	 * @param compute_r
	 * @return false if the two particles are not bonded neighbours, true otherwise
	 */
	bool _check_bonded_neighbour(BaseParticle **p, BaseParticle **q, bool compute_r);

public:
	enum {
		BACKBONE = 0,
		BONDED_EXCLUDED_VOLUME = 1,
		STACKING = 2,
		NONBONDED_EXCLUDED_VOLUME = 3,
		HYDROGEN_BONDING = 4,
		CROSS_STACKING = 5,
		COAXIAL_STACKING = 6
	};
	DNAInteraction();
	virtual ~DNAInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	bool has_custom_stress_tensor() const override {
		return true;
	}

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	// model constants
	number F1_EPS[2][5][5];
	number F1_SHIFT[2][5][5];
	number F1_A[2];
	number F1_RC[2];
	number F1_R0[2];
	number F1_BLOW[2];
	number F1_BHIGH[2];
	number F1_RLOW[2];
	number F1_RHIGH[2];
	number F1_RCLOW[2];
	number F1_RCHIGH[2];

	number F2_K[2];
	number F2_RC[2];
	number F2_R0[2];
	number F2_BLOW[2];
	number F2_RLOW[2];
	number F2_RCLOW[2];
	number F2_BHIGH[2];
	number F2_RCHIGH[2];
	number F2_RHIGH[2];

	number F4_THETA_A[13];
	number F4_THETA_B[13];
	number F4_THETA_T0[13];
	number F4_THETA_TS[13];
	number F4_THETA_TC[13];

	number F5_PHI_A[4];
	number F5_PHI_B[4];
	number F5_PHI_XC[4];
	number F5_PHI_XS[4];
};


number DNAInteraction::_repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
			if(update_forces) force = -r * (2 * EXCL_EPS * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}

#endif /* DNA_INTERACTION_H */

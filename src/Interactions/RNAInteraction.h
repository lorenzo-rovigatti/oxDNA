/*
 * BaseInteraction is a class that is meant "contain" all the interaction
 * types
 */


/**
 * @brief Handles interactions between RNA nucleotides.
 *
 * It supports both the average-base or sequence-dependent parametrization 
 * average-base is the default interaction.
 * 
 * If you want to use the RNA model, you need to set
 *  interaction_type = RNA in input file
 * Input options:
 *
 * @verbatim
[use_average_seq = <boolean> (defaults to yes)]
[seq_dep_file = <string> (sets the location of the files with sequence-dependent parameters)]
[external_model = <string> (overrides default constants for the model, set in rna_model.h), by values specified by this option)]
@endverbatim
 */

#ifndef RNA_INTERACTION_H
#define RNA_INTERACTION_H

#include "BaseInteraction.h"
#include "rna_model.h"

class RNAInteraction : virtual public BaseInteraction {
protected:
	virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	inline number _repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);

protected:
	bool _average;

	/// variables used when max_backbone_force = true
	bool _use_mbf;
	number _mbf_xmax;
	number _mbf_fmax;
	number _mbf_finf;

	char _seq_filename[512];

	Model *model;
    number _T;
	number _cross_seq_dep_K[5][5];

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

	/**
	 * @brief Check the relation between p and q. Used by the bonded interaction terms.
	 *
	 * Check whether q is the 3' neighbour of p. If this is not the case, it changes identities
	 * between q and p if p is the 3' neighbour of q. It return false if the two particles are
	 * not bonded neighbours, true otherwise.
	 * @param p
	 * @param q
	 * @param r
	 *
	 * @return false if the two particles are not bonded neighbours, true otherwise
	 */
	bool _check_bonded_neighbour(BaseParticle **p, BaseParticle **q, bool compute_r);

	bool _are_bonded(BaseParticle *p, BaseParticle *q) { return (p->n3 == q || p->n5 == q); }

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
	RNAInteraction();
	virtual ~RNAInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();
	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	Model *get_model() {return model;}

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);

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

	number F4_THETA_A[20];
	number F4_THETA_B[20];
	number F4_THETA_T0[20];
	number F4_THETA_TS[20];
	number F4_THETA_TC[20];

	number F5_PHI_A[4];
	number F5_PHI_B[4];
	number F5_PHI_XC[4];
	number F5_PHI_XS[4];
};

#endif /*RNA_INTERACTION_H*/

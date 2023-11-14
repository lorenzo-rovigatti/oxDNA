
#ifndef DRH_INTERACTION_H
#define DRH_INTERACTION_H

#include "DNA2Interaction.h"
#include "RNAInteraction2.h"
#include "rna_model.h"
#include "drh_model.h"
#include "../Particles/RNANucleotide.h"
#include "../Particles/DNANucleotide.h"

/**
 * Eryk 2022/23
 * 
 * 
 * Handles nucleic acid interactions in DNA/RNA hybrid (DRH) systems. As much as possible
 * is inherited from DNA2/RNA2, and new hybrid inter-strand interactions are defined here.
 * 
 * To use, set interaction_type = NA.
 *
 * 
 **/

class DRHInteraction : virtual public DNA2Interaction, virtual public RNA2Interaction {
protected:
	std::string _seq_filename;
	std::vector<char> _acid_type;

	number _f4_DRH(number t, int type);
	number _f4D_DRH(number t, int type);
	number _f4Dsin_DRH(number t, int type);
	number _f5_DRH(number f, int type);
	number _f5D_DRH(number f, int type);

	number _fakef4_DRH(number t, void * par);
	number _fakef4D_DRH(number t, void * par);
	number _fakef4_cxst_t1_DRH(number t, void * par);
	number _fakef4D_cxst_t1_DRH(number t, void * par);

	
	number _f1_DRH(number r, int type, int n3, int n5, bool is_n3_DNA, bool is_n5_DNA);
	number _f1D_DRH(number r, int type, int n3, int n5, bool is_n3_DNA, bool is_n5_DNA);
	number _f2_DRH(number r, int type);
	number _f2D_DRH(number r, int type);

	inline number _repulsive_lj_DRH(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);

	// DRH-only interactions
	virtual number _hydrogen_bonding_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _nonbonded_excluded_volume_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking_DRH(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	// DNA/RNA/DRH interactions
	virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);



	//These are inherited from DNA2Interaction:

	/**
	 * @brief Custom function that returns f4. This was added to add the possibility to avoid the use of meshes in classes that inherit from this.
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction (which mesh to use)
	 */
	//virtual number _custom_f4_DRH (number cost, int i) { return this->_mesh_f4[i].query(cost); }
	virtual number _custom_f4_DRH (number cost, int i) { return DNA2Interaction::_custom_f4(cost, i); }
	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	//virtual number _custom_f4D_DRH (number cost, int i) { return this->_mesh_f4[i].query_derivative(cost); }
	virtual number _custom_f4D_DRH (number cost, int i) { return DNA2Interaction::_custom_f4D(cost, i); }


	//CHANGES TO .cpp: changed the above two funcs to dna-inherited versions


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
	bool _are_bonded(BaseParticle *p, BaseParticle *q) { return (p->n3 == q || p->n5 == q); }

	bool _is_DNA(BaseParticle *p);
	int _interaction_type(BaseParticle *p, BaseParticle *q);



public:
	enum {
		BACKBONE = 0,
		BONDED_EXCLUDED_VOLUME = 1,
		STACKING = 2,
		NONBONDED_EXCLUDED_VOLUME = 3,
		HYDROGEN_BONDING = 4,
		CROSS_STACKING = 5,
		COAXIAL_STACKING = 6,
		DEBYE_HUCKEL = 7
	};


	DRHInteraction();
	virtual ~DRHInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();  //apparently renaming many or all of the below methods causes errors
	virtual void allocate_particles(std::vector<BaseParticle *> &particles);  

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	
	// DRH model constants
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


#endif 


#ifndef DNA_WITH_RNA_INTERACTION_H
#define DNA_WITH_RNA_INTERACTION_H

#include "BaseInteraction.h"
#include "rna_model.h"
#include "hybrid_model.h"
#include "../Particles/RNANucleotide.h"

/**
 * @brief Handles interactions between DNA and RNA nucleotides.
 *
 * 
 */

class DNAwithRNAInteraction : public BaseInteraction {
protected:
	std::string _nucleotide_types;

	bool _average_DNA;
	bool _average_RNA;
	bool _average_hybrid;

	std::string _seq_filename_DNA;
	char _seq_filename_RNA[512];
	char _seq_filename_hybrid[512];

	number _T;
	number _cross_seq_dep_K_RNA[5][5];
	number _cross_K_multiplier_hybrid;
	number _hb_multiplier_DNA;
	bool _grooving_DNA;

	/// true by default; set this to false if you want the code to not die when bonded backbones are found to be outside the acceptable FENE range
	bool _allow_broken_fene_DNA;

	/// here we store the r0, which is different in oxDNA and oxDNA2, and is set in the constructor
	number _fene_r0_DNA;

	/// variables used when max_backbone_force = true
	bool _use_mbf_DNA;
	number _mbf_xmax_DNA;
	number _mbf_fmax_DNA;
	number _mbf_finf_DNA;

	bool _use_mbf_RNA;
	number _mbf_xmax_RNA;
	number _mbf_fmax_RNA;
	number _mbf_finf_RNA;


	int MESH_F4_POINTS_DNA[13];
	Mesh _mesh_f4_DNA[13];
	int MESH_F4_POINTS_RNA[13];
	Mesh _mesh_f4_RNA[13];
	int MESH_F4_POINTS_hybrid[13];
	Mesh _mesh_f4_hybrid[13];

	//(for RNA)
	Model *model;

	
	virtual void _on_T_update(); 


	number _f1_DNA(number r, int type, int n3, int n5);
	number _f1D_DNA(number r, int type, int n3, int n5);
	number _f2_DNA(number r, int type);
	number _f2D_DNA(number r, int type);
	number _f4_DNA(number t, int type);
	number _f4D_DNA(number t, int type);
	number _f4Dsin_DNA(number t, int type);
	number _f5_DNA(number f, int type);
	number _f5D_DNA(number f, int type);

	number _fakef4_DNA(number t, void * par);
	number _fakef4D_DNA(number t, void * par);
	number _fakef4_cxst_t1_DNA(number t, void * par);
	number _fakef4D_cxst_t1_DNA(number t, void * par);

	number _f1_RNA(number r, int type, int n3, int n5);
	number _f1D_RNA(number r, int type, int n3, int n5);
	number _f2_RNA(number r, int type);
	number _f2D_RNA(number r, int type);
	number _f4_RNA(number t, int type);
	number _f4D_RNA(number t, int type);
	number _f4Dsin_RNA(number t, int type);
	number _f5_RNA(number f, int type);
	number _f5D_RNA(number f, int type);

	number _fakef4_RNA(number t, void * par);
	number _fakef4D_RNA(number t, void * par);
	number _fakef4_cxst_t1_RNA(number t, void * par);
	number _fakef4D_cxst_t1_RNA(number t, void * par);

	//for hybrid h-bonding
	number _f1_hybrid(number r, int type, int n3, int n5, bool is_n3_DNA, bool is_n5_DNA);
	number _f1D_hybrid(number r, int type, int n3, int n5, bool is_n3_DNA, bool is_n5_DNA);

	//for hybrid cross-stacking
	number _f2_hybrid(number r, int type);
	number _f2D_hybrid(number r, int type);


	virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);



	virtual number _backbone_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	virtual number _nonbonded_excluded_volume_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);


	virtual number _backbone_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _bonded_excluded_volume_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	virtual number _nonbonded_excluded_volume_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking_RNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	inline number _repulsive_lj_RNA(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);
	inline number _repulsive_lj_DNA(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);

	virtual number _hydrogen_bonding_hybrid(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _nonbonded_excluded_volume_hybrid(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking_hybrid(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _coaxial_stacking_hybrid(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	/**
	 * @brief Custom function that returns f4. This was added to add the possibility to avoid the use of meshes in classes that inherit from this.
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4_DNA (number cost, int i) { return this->_mesh_f4_DNA[i].query(cost); }

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4D_DNA (number cost, int i) { return this->_mesh_f4_DNA[i].query_derivative(cost); }


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
		COAXIAL_STACKING = 6
	};


	DNAwithRNAInteraction();
	virtual ~DNAwithRNAInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();  //apparently renaming many or all of the below methods causes errors
	virtual void allocate_particles(std::vector<BaseParticle *> &particles);  

	Model *get_model() {return model;}

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	// DNA model constants
	number F1_EPS_DNA[2][5][5];
	number F1_SHIFT_DNA[2][5][5];
	number F1_A_DNA[2];
	number F1_RC_DNA[2];
	number F1_R0_DNA[2];
	number F1_BLOW_DNA[2];
	number F1_BHIGH_DNA[2];
	number F1_RLOW_DNA[2];
	number F1_RHIGH_DNA[2];
	number F1_RCLOW_DNA[2];
	number F1_RCHIGH_DNA[2];

	number F2_K_DNA[2];
	number F2_RC_DNA[2];
	number F2_R0_DNA[2];
	number F2_BLOW_DNA[2];
	number F2_RLOW_DNA[2];
	number F2_RCLOW_DNA[2];
	number F2_BHIGH_DNA[2];
	number F2_RCHIGH_DNA[2];
	number F2_RHIGH_DNA[2];

	number F4_THETA_A_DNA[13];
	number F4_THETA_B_DNA[13];
	number F4_THETA_T0_DNA[13];
	number F4_THETA_TS_DNA[13];
	number F4_THETA_TC_DNA[13];

	number F5_PHI_A_DNA[4];
	number F5_PHI_B_DNA[4];
	number F5_PHI_XC_DNA[4];
	number F5_PHI_XS_DNA[4];



	// RNA model constants
	number F1_EPS_RNA[2][5][5];
	number F1_SHIFT_RNA[2][5][5];
	number F1_A_RNA[2];
	number F1_RC_RNA[2];
	number F1_R0_RNA[2];
	number F1_BLOW_RNA[2];
	number F1_BHIGH_RNA[2];
	number F1_RLOW_RNA[2];
	number F1_RHIGH_RNA[2];
	number F1_RCLOW_RNA[2];
	number F1_RCHIGH_RNA[2];

	number F2_K_RNA[2];
	number F2_RC_RNA[2];
	number F2_R0_RNA[2];
	number F2_BLOW_RNA[2];
	number F2_RLOW_RNA[2];
	number F2_RCLOW_RNA[2];
	number F2_BHIGH_RNA[2];
	number F2_RCHIGH_RNA[2];
	number F2_RHIGH_RNA[2];

	number F4_THETA_A_RNA[20];
	number F4_THETA_B_RNA[20];
	number F4_THETA_T0_RNA[20];
	number F4_THETA_TS_RNA[20];
	number F4_THETA_TC_RNA[20];

	number F5_PHI_A_RNA[4];
	number F5_PHI_B_RNA[4];
	number F5_PHI_XC_RNA[4];
	number F5_PHI_XS_RNA[4];


	// Hybrid constants
	number F1_EPS_hybrid[2][5][5];
	number F1_SHIFT_hybrid[2][5][5];
	number F1_A_hybrid[2];
	number F1_RC_hybrid[2];
	number F1_R0_hybrid[2];
	number F1_BLOW_hybrid[2];
	number F1_BHIGH_hybrid[2];
	number F1_RLOW_hybrid[2];
	number F1_RHIGH_hybrid[2];
	number F1_RCLOW_hybrid[2];
	number F1_RCHIGH_hybrid[2];

	number F2_K_hybrid[2];
	number F2_RC_hybrid[2];
	number F2_R0_hybrid[2];
	number F2_BLOW_hybrid[2];
	number F2_RLOW_hybrid[2];
	number F2_RCLOW_hybrid[2];
	number F2_BHIGH_hybrid[2];
	number F2_RCHIGH_hybrid[2];
	number F2_RHIGH_hybrid[2];

	number F4_THETA_A_hybrid[13];
	number F4_THETA_B_hybrid[13];
	number F4_THETA_T0_hybrid[13];
	number F4_THETA_TS_hybrid[13];
	number F4_THETA_TC_hybrid[13];

	number F5_PHI_A_hybrid[4];
	number F5_PHI_B_hybrid[4];
	number F5_PHI_XC_hybrid[4];
	number F5_PHI_XS_hybrid[4];
};


#endif /* DNA_INTERACTION_H */

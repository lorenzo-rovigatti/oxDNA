/*
 * PatchyShapeInteraction.h
 *
 *  Created on: 14/jul/2018
 *      Author: petr
 */

#ifndef PATCHYSHAPEINTERACTION_H_
#define PATCHYSHAPEINTERACTION_H_


#define PLPATCH_VM1 0
#define PLPATCH_VM3 1

#define PATCHY_POWER 200
#define PATCHY_CUTOFF 0.18f

#define PLEXCL_S   1.0f
#define PLEXCL_R   0.9053f
#define PLEXCL_B   677.505671539f
#define PLEXCL_RC  0.99888f
#define PLEXCL_EPS 2.0f

//constants based on villars paper
#define PLEXCL_VS   1.0f
#define PLEXCL_VR   1.15f
#define PLEXCL_VB   (-0.3797564833966851)
#define PLEXCL_VRC  (2.757802724660435)
#define PLEXCL_VEPS 2.0f

//48-96 lj potnetial
#define PLEXCL_V48S   1.0
#define PLEXCL_V48R   1.15f
#define PLEXCL_V48B   (-2.11835)
#define PLEXCL_V48RC  1.19798
#define PLEXCL_V48EPS 2.0f

//there are two types of modulation
#define PLEXCL_NARROW_N 2

#include "../../../../src/Interactions/BaseInteraction.h"
#include "../Particles/PatchyShapeParticle.h"
#include "../../../../src/Observables/BaseObservable.h"
/**
 * @brief Manages the interaction between simple patchy particles, each can have multiple patches of different colors and rigid body of different shapes (
 *
 * This interaction is selected with
 * interaction_type = plpatchy
 *
 * @verbatim
 *
particle_types_N = <int> (number of particle types)
patch_types_N  = <int> (number of patches type
patchy_file = <string> (name of the file with patch types)
particle_file = <string> (name of the file with particle types)
use_torsion = <bool> (use or not use torsion constraints for binding, default to 1)
same_type_bonding = <bool> (two particles of the same type can bind, default to 1)
narrow_type = <int> (for lj48 like interaction, sets the type of narrowness of the bonds)
interaction_tensor = <bool> (false by default; if true, possible interactions are loaded from a file specified by the option below)
interaction_tensor_file = <string> (filename of the interaction tensor file; interactions specified in a following way described below)
same_type_bonding = <bool> (particles of the same type can bond)
no_multipatch = <bool> (if set to 1, the code does not allow 1 patch binding to more than 1 other patch, and uses the lock patch; only works if used with MC2 and special move MCPatchyShapeMove

 input config files about types of patches will be in a separate file of the following form:
	 *   patch_0 = {
	 *    id = some integer, preferably the same as in the case of patch name
	 *    color = some integer (only patches of the same color interact)
	 *    strength = 1.0 (strength of the patch-patch interaction; The same color patches need to have the same strength
	 *    position = x,y,z   position of the patch wrt to the c.o.m of the particle
	 *    a1 = x,y,z      a1 and a2 specify the orientation of the patch
	 *    a2 = x,y,z
	 *   }
	 *
	 *   particle_0 = {
	 *    id = some_integer
	 *    N_patches = integer
	 *    patches = 1,2,4,5,6,8,7,8
	 *   }
	 *
if tensor file is used, the format of allowed interactions is:
particle_type_id1 particle_type_id2  patch_1  patch_2

@endverbatim
 */

enum TorsionVersion {
    NO_TORSION, // no torsional modulation
    PARTIAL_TORSION, // only a1 vector angle used (new addition)
    FULL_TORSION // a1 and a3 vectors used
};

class PatchyShapeInteraction: public BaseInteraction {
protected:
	/// Number of patches per particle

//    int _N_patch_types; //the total number of different patch types; each patchy particle can have multiple patches, of different types
//    int _N_particle_types; //number of different particle types
    int N_patches; //total number of patches in the simulation

	/// Repulsive interaction energy at the cut-off
	number _E_cut;

	/// Patch-patch interaction energy at the cut-off
	number _patch_E_cut;

	/// Width of the patch, defaults to 0.12
	number _patch_alpha;

	/// _patch_alpha^10
	number _patch_pow_alpha;

	int _shape;
	int _patch_type;

	number _kf_delta, _kf_cosmax;


	int _narrow_type; //the type of angular narrowness used by the patches
    TorsionVersion _use_torsion;

    bool _same_type_bonding; //are particles of the same type allowed to bond, by default yes

    bool _interaction_tensor;

    int *_interaction_table_types; //a 2d matrix of _N_particle_types * _N_particle_types that checks if 2 particles can interact (filled with 0 and 1)
    //int *_interaction_table_types;
    int *_interaction_patch_types; //a 2d matrix of patches (each patch has a unique id) that checks if patches can interact; is only used if interaction tensor is true

    bool _no_multipatch;

    number _lock_cutoff;

	std::vector<Patch> _patch_types;
	std::vector<PatchyShapeParticle> _particle_types;

	number _sphere_radius;

	//model constants
	number PLPATCHY_THETA_TC[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_T0[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_TS[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_A[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_B[PLEXCL_NARROW_N];

	//icosahedron constants
	int *_close_vertexes; // 5 close vertexes for each vertex
	number _tworinscribed ; // twice the radius of inscribed sphere


public:
	virtual bool _bonding_allowed(PatchyShapeParticle  *p, PatchyShapeParticle  *q, int pi, int pj );
	bool _patches_compatible(PatchyShapeParticle  *p, PatchyShapeParticle  *q, int pi, int pj );

    int num_patch_types() const {return _patch_types.size();};
    int num_particle_types() const {return _particle_types.size();};

	number _V_mod(int type, number cosr1);
	number _V_modD(int type, number cosr1);
	number _V_modDsin(int type, number cosr1);
	number _exc_vol_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);

	number _repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);
	number _repulsive_lj_n(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc,int n,bool update_forces);

	number _exc_LJ_vol_interaction_sphere(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
//	number _exc_vol_hard_icosahedron(BaseParticle *p, BaseParticle *q, LR_vector *r,bool update_forces);
	number _exc_vol_hs (BaseParticle *ap, BaseParticle *aq, bool compute_r, bool update_forces);
	//inline number _exc_quadratic_vol_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);

	inline number patchy_interaction_full_torsion(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	inline number patchy_interaction_partial_torsion(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	inline number patchy_interaction_no_torsion(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	inline number _patchy_interaction_kf(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	virtual bool multipatch_allowed() {return  ! this->_no_multipatch;}

	//true if particles of of this type can interact // todo: use a better container
    // todo reenable interaction tensors
//    bool can_interact(BaseParticle *p, BaseParticle *q) {return (bool) this->_interaction_table_types[p->type*  num_particle_types() + q->type]; }

	//virtual int check_valence(ConfigInfo &conf_info) {return 0;} //scans all interacting particles if one patch is bond to more particles, it breaks all bonds but 1;
	virtual number just_two_patch_interaction(PatchyShapeParticle *p, PatchyShapeParticle *q, int pi,int  qi,LR_vector *r);

	number get_patch_cutoff_energy() const {return this->_lock_cutoff;}

	Patch _process_patch_type(std::string input_string); //this function processes patch type from the input file
	PatchyShapeParticle _process_particle_type(std::string input_string);

    void _load_patchy_particle_files(std::string& patchy_file, std::string& particle_file);

    void _load_interaction_tensor(std::string &tensor_file); //only used if tensor file provided

    void _init_icosahedron();

    void _init_patchy_locks(std::shared_ptr<ConfigInfo> _Info = NULL);

    void check_patchy_locks(std::shared_ptr<ConfigInfo> _Info = NULL);

public:
	enum {
		PATCHY = 0,
		EXCVOL = 1,
	};

	enum {
		SPHERE_SHAPE = 0,
		ICOSAHEDRON_SHAPE = 1,
		HS_SHAPE = 2,
	};

	enum {
		POINT_PATCHES = 0,
		KF_PATCHES = 1,
	};

	PatchyShapeInteraction();
	virtual ~PatchyShapeInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	number get_alpha() const { return _patch_alpha; }


//	void check_loaded_particles(); //needed for debugging

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r=false, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r= false, bool update_forces= false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r=false, bool update_forces=false);
//	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r=NULL, bool update_forces=false) {
//		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
//	}

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle*> &particles);

	//virtual void generate_random_configuration(BaseParticle **particles, int N, number box_side);
};

#endif /* PATCHYINTERACTION_H_ */

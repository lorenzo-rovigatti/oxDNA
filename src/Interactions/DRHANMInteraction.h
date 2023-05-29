/*
 * DRHANMInteraction.h
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 */

#ifndef DRHANM_INTERACTION_H_
#define DRHANM_INTERACTION_H_

#include "DRHInteraction.h"
#include "BaseInteraction.h"

#include "../Particles/ACParticle.h"
#include "../Particles/DNANucleotide.h"
#include "../Particles/RNANucleotide.h"
#include "rna_model.h"


class DRHANMInteraction: public DRHInteraction {

protected:
	//since the new topology doesn't support proteins (yet)
	std::string _nucleotide_types;
	
	int ndna;//How many particles of DNA type: Used in allocate_particles
	int npro;//How many bonds b/t different particle types
	int ndnas;//Number of Strands that are dna
	int _firststrand; //+ for dna, - for protein

	number _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    number _pro_base_sigma,_pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;
	number _pro_sigma, _pro_rstar, _pro_b, _pro_rcut;

	std::map<std::pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
	std::map<std::pair<int, int>, std::pair<char, double> > _potential;

	//map<pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
	//map<pair<int, int>, pair<char, double> > _potential;


public:
	enum {
			SPRING = 8,
			PRO_EXC_VOL = 9,
			PRO_DNA_EXC_VOL = 10
			//Assigned 8 9 and 10 so it won't overwrite the already existing DNA function pointers in the _int_map
	};

    char _parameterfile[500];

	DRHANMInteraction();
	virtual ~DRHANMInteraction();
	virtual void get_settings(input_file &inp); //done

	virtual void allocate_particles(std::vector<BaseParticle*> &particles); 
	virtual void read_topology(int *N_strands, std::vector<BaseParticle*> &particles); 
		//^these two may need changing

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r=true, bool update_forces=false); //done
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false); //done
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false); //done

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles); //done
	virtual void init(); //done


	number _protein_na_exc_volume(BaseParticle *p,BaseParticle *q, LR_vector *r, bool update_forces); //done
	number _protein_na_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness); 
	number _protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces); //done      ^this one done also
	number _protein_exc_volume(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces); //done
	virtual number _protein_spring(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces); //done


    virtual number _na_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _na_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done

    //for the newer code the args (in declarations) for interactions are: 'BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces'
    //when they are called the args should be 'p, q, false, update_forces'
};

#endif 

/*
 * DNANMInteraction.h
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 */

#ifndef DNANM_INTERACTION_H_
#define DNANM_INTERACTION_H_

#include "DNA2Interaction.h"
#include "BaseInteraction.h"

#include "../Particles/ACParticle.h"
#include "../Particles/DNANucleotide.h"


class DNANMInteraction: public DNA2Interaction {

protected:
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

	DNANMInteraction();
	virtual ~DNANMInteraction();
	virtual void get_settings(input_file &inp); //done

	virtual void allocate_particles(std::vector<BaseParticle*> &particles); //done
	virtual void read_topology(int *N_strands, std::vector<BaseParticle*> &particles); //done


	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r=true, bool update_forces=false); 
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false); 
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false); 

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
	virtual void init();


	number _protein_dna_exc_volume(BaseParticle *p,BaseParticle *q, LR_vector *r, bool update_forces); //original copied
	number _protein_dna_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness); //original copied
	number _protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces);
	number _protein_exc_volume(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);
	virtual number _protein_spring(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);


    virtual number _dna_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done
    virtual number _dna_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces); //done

    //for the newer code the args (in definitions) for interactions are: 'BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces'
    //when they are called the args should be 'p, q, false, update_forces'
};

#endif /* DNANM_INTERACTION_H_ */

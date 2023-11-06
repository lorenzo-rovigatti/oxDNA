// Patchy Helix Bundles

#ifndef PHBInteraction_H_
#define PHBInteraction_H_
#include "BaseInteraction.h"
#include "../Particles/PHBParticle.h"
#include <sstream>


class PHBInteraction: public BaseInteraction{
public:
    PHBInteraction();
    virtual ~PHBInteraction();
    // Variables d
    int totPar,strands,totHelix,totPatchy,i; //topology
    bool _allow_broken_fene = false,_prefer_harmonic_over_fene = false; //special implementation for later
    std::string temp;
    // Necessary interaction
    virtual void get_settings(input_file &inp);
	virtual void init();
	virtual void allocate_particles(std::vector<BaseParticle *> &particles); //Add particle to the system
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles); // Read the top file
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles); // Check all the input file are correct.

	//Interaction that are updated repeatedly
	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false); //Check bonded or non-bonded
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false); // Bonded particle interaction
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false); //Non-bonded particle interaction
    
    //Common Interaction
    virtual number _repulsive_lj2(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces); //Excluded volume with pre-factor
    virtual number maxRadius(std::vector<PHBParticle*> &particles);//returns the maximum radius
    //Helix Interactions
    virtual number spring(BaseParticle *p, BaseParticle *q, bool compute_r=true, bool update_forces=false); //Calculate spring interaction



protected:
};

























#endif
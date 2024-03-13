// This is 2nd itteration of PHB interaction

#ifndef PHBINTERACTION2_H
#define PHBINTERACTION2_H

#include "BaseInteraction.h"
#include "../Particles/CCGParticle.h"

class PHBInteraction2: public BaseInteraction {
public:
    // Temporary variables
    int i,j;
    // Patchy variables



    PHBInteraction2();
    virtual ~PHBInteraction2();

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


}

#endif // PHBINTERACTION2_H
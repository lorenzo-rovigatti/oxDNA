// CCGInteraction Header
// Subhajit
/*
** btype = color of the particles
*/


#ifndef CCGInteraction_H_
#define CCGInteraction_H_

#include "BaseInteraction.h"
#include "../Particles/PatchyParticle.h"
#include <fstream>
#include <sstream>
class CCGInteraction: public BaseInteraction {
protected:
public:
	int version,totPar,strands,ccg,ccg0,noSpring,noColor ;// Header Information
	int currentVersion = 1; // Don't forget to update the version number
	int particleType,color,neighbour,bfactor; //Body parameters particleType,particleName,...
	double Bfactor;
	// enum {
	// 	PATCHY = 4
	// };
	int i,j;
	std::string temp;
	bool connection; // connection shifts between adding spring neighbours and Bfactor during reading of the topology file
	CCGInteraction();
	virtual ~CCGInteraction();

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

	//My interactions
	virtual number spring(BaseParticle *p, BaseParticle *q, bool compute_r=true, bool update_forces=false); //Calculate spring interaction
	virtual number exc_vol(BaseParticle *p, BaseParticle *q, bool compute_r=true,bool update_forces=false); //Calculate excluded volume interaction
};

#endif /* CCGInteraction_H_ */
// CCGInteraction Header
// Subhajit
/*
** btype = color of the particles
*/


#ifndef CCGInteraction_H_
#define CCGInteraction_H_

#include "BaseInteraction.h"
#include "../Particles/PatchyParticle.h"
#include "../Utilities/parse_input/parse_input.h"
#include <fstream>
#include <sstream>
class CCGInteraction: public BaseInteraction {
protected:
public:
	int version,totPar,strands,ccg,ccg0,noSpring,noColor ;// Header Information
	int currentVersion = 1; // Don't forget to update the version number
	int particleType,color,neighbour,bfactor; //Body parameters particleType,particleName,...
	double patchyRcut=1.2,patchyAlpha=0.12,patchyRadius=0,patchyCutoff=0;//color parameters
	double Bfactor,strength =1.0,rmod;
	LR_vector r; //temporary parameters.
	// enum {
	// 	PATCHY = 4
	// };
	int i,j;
	std::string temp;
	bool connection,bcall; // connection shifts between adding spring neighbours and Bfactor during reading of the topology file
	const double sigma=1.0f,rstar=0.9053f,b=677.505671539f,rc=0.99888f,epsilon=2.0f;
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

	//Color interactions
	virtual bool color_compatibility(BaseParticle *p, BaseParticle *q); //check wether two particle will interact or not
	virtual number patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r=true,bool update_forces=false);
};

#endif /* CCGInteraction_H_ */
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

    virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

#endif /* CCGInteraction_H_ */
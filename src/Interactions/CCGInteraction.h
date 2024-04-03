// CCGInteraction Header
// Subhajit
/*
** btype = color of the particles
*/


#ifndef CCGInteraction_H_
#define CCGInteraction_H_

// Make sure that PHB is also set to greater or same value
#ifndef MAXparticles
#define MAXparticles 2000
#endif
#ifndef MAXneighbour
#define MAXneighbour 13
#endif


#include "BaseInteraction.h"
#include "../Particles/PatchyParticle.h"
#include "../Utilities/parse_input/parse_input.h"
#include <fstream>
#include <sstream>
#include <cmath>
class CCGInteraction: public BaseInteraction {
protected:
public:
	enum {
		CCG = 0
	};
	int version,totPar,strands,ccg,ccg0,noSpring,noColor ;// Header Information
	int currentVersion = 1; // Don't forget to update the version number
	int particleType,color,neighbour,bfactor; //Body parameters particleType,particleName,...
	double patchyRcut=1.2,patchyAlpha=0.12,patchyRadius=0,patchyCutoff=0.2,patchyEcutoff=0,patchyB=667.505671539f;//color parameters // cut off should be 0.0324
	double Bfactor,strength =1,rmod,rnorm;
	double damp=1.0;
	// LR_vector r; //temporary parameters.
	enum {
		SPRING =1,
		EXEVOL=2,
		PATCHY=3,
		EXEVOLN=4
	};
	int i,j;
	std::string temp;
	// bool connection,bcall; // connection shifts between adding spring neighbours and Bfactor during reading of the topology file
	
	double patchySigma=1.0f,patchyRstar=0.9053f,patchyRc=0.99998,patchyEpsilon=2.0f,patchyLockCutOff=0,patchyInteractionDistanceCutoff=0;

	//GPU matrix
	int CPUconnections[MAXparticles][MAXneighbour+1];//first intger will state number of connections
	double CPUro[MAXparticles][MAXneighbour+1];//rarius of the spring
	double CPUk[MAXparticles][MAXneighbour+1];//spring constant
	int CPUtopology[MAXparticles][4];// strand, strength, iC,radius
	// double patches[GPUmaxiP][5];// color,strength,x,y,z


	CCGInteraction();
	virtual ~CCGInteraction();

	number ccg_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false){
		OX_DEBUG("Bonded interaction is called");
		number energy =0.f;
		
		energy+=25;
		return energy;
	};
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
	virtual number exc_vol_bonded(BaseParticle *p, BaseParticle *q, bool compute_r=true,bool update_forces=false); //Calculate excluded volume interaction
	virtual number exc_vol_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r=true,bool update_forces=false); //Calculate excluded volume interaction
	// virtual number exc_vol_centralParticle(BaseParticle *p, BaseParticle *q, bool compute_r=true, bool update_forces=false); // Special calculation for base particle

	//Color interactions
	virtual bool color_compatibility(BaseParticle *p, BaseParticle *q); //check wether two particle will interact or not
	virtual number patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r=true,bool update_forces=false);

	//Debug function
	virtual number debug(BaseParticle  *p, BaseParticle*q, bool compute_r,bool update_forces);
	virtual number _repulsive_lj(const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces);
};

#endif /* CCGInteraction_H_ */
// This is 2nd itteration of PHB interaction

#ifndef PSPInteraction_H
#define PSPInteraction_H

//Model maximum capacity
#define PSPmaxParticles 2000
#define PSPmaxNeighbour 20
#define PSPmaxParticleColor 20
#define PSPmaxPatchColor 20
#define PSPmaxPatchOnParticle 6 // For PHB interaction only

#include "BaseInteraction.h"
#include "../Particles/CCGParticle.h"
#include <sstream>
#include "omp.h" // for openmp

class PSPInteraction: public BaseInteraction {
public:
    // Temporary variables
    int i,j;
	number rmod, rnorm;
	std::string temp,line;


	//Header variables
	int totPar,strads,totParticleColor,totPatchColor; //topology
	bool harmonics=true;

	// Patchy Parameters
	double patchyRcut=1.2,patchyAlpha=0.12,patchyRadius=0,patchyCutoff=0.2,patchyEcutoff=0,patchyB=667.505671539f;

	//Topology variables
	int connections[PSPmaxParticles][PSPmaxNeighbour]; // 5 0 1 2 3 4 where 5 is the number of neighbours and 0 1 2 3 4 are the neighbours
	float r0[PSPmaxParticles][PSPmaxNeighbour]; // Equilibrium radius of the spring, if starts with 0 all particles have different radius, if 1 all particles have same radius
	float k0[PSPmaxParticles][PSPmaxNeighbour]; // Spring constant, same as above
	float particleTopology[PSPmaxParticles][3]; // Strand, particleColor, radius
	float patches[PSPmaxPatchColor][5]; // color, strength, x, y, z // for PSP interaction patch color is useless and should be -1
	int particlePatches[PSPmaxParticleColor][PSPmaxPatchColor]; // Number of patches, patch1, patch2, patch3, patch4, patch5, patch6
    // Patchy variables



    PSPInteraction();
    virtual ~PSPInteraction();

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


	//Custom Interactions
	virtual number spring(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false); //Spring interaction
	virtual number exeVolInt(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false); //Excluded volume interaction
	virtual number linearRepulsion(number patchyEpsilon, LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces); //Linear repulsion
	virtual number cubicRepulsion(number patchyEpsilon, LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces); //Cubic repulsion

	//Checkers
	int bonded(int i, int j);


};

#endif // PSPInteraction_H
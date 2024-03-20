// This is 2nd itteration of PHB interaction

#ifndef PSPInteraction_H
#define PSPInteraction_H

//Model maximum capacity
#define PSPmaxParticles 20
#define PSPmaxNeighbour 20
#define PSPmaxParticleColor 20
#define PSPmaxPatchColor 20
#define PSPmaxPatchOnParticle 6 // For PHB interaction only

#include "BaseInteraction.h"
#include "../Particles/CCGParticle.h"
#include <sstream>
#include <fstream>
#include <cmath>
#include "omp.h" // for openmp

class PSPInteraction: public BaseInteraction {
private:
	//Patch variables
	int patchType[PSPmaxPatchColor];
	// void coordinateConverter();
public:
    // Temporary variables
    int i,j;
	number rmod, rnorm;
	std::string temp,line;

	// Interaction Modes
	bool helixBubble=true,harmonics=true;

	//Header variables
	int totPar,strads,totParticleColor,totPatchColor; //topology

	// Patchy Parameters
	number patchyRcut=1.2,patchyAlpha=0.12,patchyRadius=0,patchyCutoff=0.2,patchyEcutoff=0,patchyB=667.505671539f,patchySigma=1.0f,patchyRstar=0.9053f,patchyRc=0.99998,patchyEpsilon=2.0f,patchyEcut=0,patchyLockCutOff=-0.1;
	number patchyRcut2=SQR(patchyRcut),patchyPowAlpha = powf(patchyAlpha, (number) 10.f);
	number tepEpsilon=1.0f,tepB=1,_xu_bending=0.952319757,_xk_bending= 1.14813301,tepFeneDelta=1.6,_ka = 100.,_kb = 21.,_kt = 29.7,_twist_a = 0.00,_twist_b = 0.95;
	number tepFeneDelta2=SQR(tepFeneDelta);

	//Topology variables
	int connections[PSPmaxParticles][PSPmaxNeighbour+1]; // 5 0 1 2 3 4 where 5 is the number of neighbours and 0 1 2 3 4 are the neighbours
	float r0[PSPmaxParticles][PSPmaxNeighbour+1]; // Equilibrium radius of the spring, if starts with 0 all particles have different radius, if 1 all particles have same radius
	float k0[PSPmaxParticles][PSPmaxNeighbour+1]; // Spring constant, same as above
	float particleTopology[PSPmaxParticles][3]; // Strand, particleColor, radius
	float patches[PSPmaxPatchColor][11]; // color, strength, x, y, z, a1x, a1y, a1z, a2x, a2y, a2z // for PSP interaction patch color is useless and should be -1
	int particlePatches[PSPmaxParticleColor][PSPmaxPatchOnParticle+1]; // Number of patches, patch1, patch2, patch3, patch4, patch5, patch6
    // Patchy variables
	bool patchLock[PSPmaxParticles][PSPmaxPatchOnParticle]; // 0 - not locked, 1 - locked


    PSPInteraction();
    virtual ~PSPInteraction();

	//Checkers
	int bonded(int i, int j);
	void setRcut(std::vector<BaseParticle *> &particles);

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
	virtual number exc_vol_nonbonded(CCGParticle *p, CCGParticle *q, bool compute_r = true,bool update_forces=false); //Excluded volume interaction
	virtual number linearRepulsion(number patchyEpsilon, LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces); //Linear repulsion
	virtual number cubicRepulsion(number patchyEpsilon, LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces); //Cubic repulsion
	virtual number fene(CCGParticle *p, CCGParticle *q, bool compute_r,bool update_forces); //Fene potential
	// Worm like potentials
	virtual number bonded_double_bending(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	virtual number bonded_alignment(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	virtual number bonded_twist(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	LR_vector rotateVectorAroundVersor(const LR_vector vector, const LR_vector versor, const number angle);

	//Patchy Interactions
	virtual number patchyInteractionSimple(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	// virtual bool patchyInteractionSimpleCheck(CCGParticle *p, CCGParticle *q);
	virtual number patchyInteractionColored(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	virtual number patchyInteraction2point(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	virtual number patchyInteraction3point(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);
	virtual number patchyInteractionBubble(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces);

	//debug
	void print2DArraytoFile(std::string filename, float *arr, int row, int col);
	void print2DArraytoFile(std::string filename, int *arr, int row, int col);


};

#endif // PSPInteraction_H


/////// Information about the topology file ///////
// Patch Types:

// 0 - Normal patch, directly the x,y,z coordinates will be used.
// 1 - Square patch, x,y,z represents the position of the patch on a unit sphere, hence the final coordiantes will be (x,y,z)*radius
// 2 - Spherical patch, where it is represented in r,theta,phi, r is not very meaningful, as it will be equal to radius of the sphere and can be any value.
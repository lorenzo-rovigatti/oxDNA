// Patchy Helix Bundles

#ifndef PHBInteraction_H_
#define PHBInteraction_H_
#include "BaseInteraction.h"
#include "../Particles/PHBParticle.h" // This is helix particles + patchy particles
#include <sstream>
#include "omp.h" // for openmp


class PHBInteraction: public BaseInteraction{
public:
    enum{
        SPRING=0,
        BONDED_DOUBLE_BENDING=1,
        BONDED_TWIST=2,
        BONDED_ALIGNMENT=3,
        PATCHY_INTERACTION_NOTORSION=4
    };
    PHBInteraction();
    virtual ~PHBInteraction();
    // Variables d
    int totPar,strands,totHelix,totPatchy,numPatches[3]={6,3,0}; //topology
    int particleType;//body parameter
    int i; //temporary parameters
    number rmod,rnorm;
    double damp=0;
    double patchySigma=1.0f,patchyRstar=0.9053f,patchyRc=0.99998,patchyEpsilon=2.0f,patchyLockCutOff=-0.1,patchyInteractionDistanceCutoff=0,patchyB=667.505671539,patchyRcut=1.2,patchyCutOff=1.8,patchyAlpha=0.12,patchyIntercept=0.0,patchySpacer=0.0;
    double patchyRcut2=SQR(patchyRcut),patchyPowAlpha = powf(patchyAlpha, (number) 10.f);
    double tepEpsilon=1.0f,tepB=1,_xu_bending=0.952319757,_xk_bending= 1.14813301,tepFeneDelta=1.6; //currently not used
    double tepFeneDelta2=SQR(tepFeneDelta);
    number _ka = 100.,_kb = 21.,_kt = 29.7; //Energies
    number 	_twist_a = 0.00,_twist_b = 0.95; //cosine

    // Helix future variables 


    bool _allow_broken_fene = false, iLoveHarmonic=true; //special implementation for later
    std::string temp;

    // Patchy variables


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
    virtual number repulsiveLJ(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces); //Excluded volume with pre-factor
    virtual number hardRepulsive(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rc, bool update_forces); //Excluded volume with power 200 potential
    virtual number exc_vol_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces);
    // virtual number maxRadius(std::vector<PHBParticle*> &particles);//returns the maximum radius

    //Helix Interactions
    virtual number fene(PHBParticle *p, PHBParticle *q, bool compute_r=true, bool update_forces=false);//Spring but with fene potential instead of harmonic potential
	virtual number spring(PHBParticle *p, PHBParticle *q, bool compute_r=true, bool update_forces=false); //Calculate spring interaction
    virtual number bonded_twist(PHBParticle *p, PHBParticle *q, bool compute_r=true, bool update_forces=false);
    virtual number bonded_double_bending(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces);
    virtual number bonded_alignment(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces);

    //Patchy Interactions
    virtual number patchy_interaction_notorsion(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces);
    virtual bool bondingAllowed(PHBParticle *p, PHBParticle *q,int pi, int qi);

    //Multiprocess Functions
    virtual void setRcut(std::vector<BaseParticle *> &particles);

    

protected:

    //Helix Interactions
    LR_vector rotateVectorAroundVersor(const LR_vector vector, const LR_vector versor, const number angle);
};

























#endif
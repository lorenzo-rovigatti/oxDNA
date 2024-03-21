// CCGParticle.h
// Subhajit

#ifndef CCGParticle_H
#define CCGParticle_H

#include "BaseParticle.h"
#include <algorithm>

struct PatchyBond {
	BaseParticle *other;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector force;
	LR_vector p_torque, q_torque;

	PatchyBond(BaseParticle *o, number my_r_p, int pp, int qp, number e) :
					other(o),
					r_p(my_r_p),
					p_patch(pp),
					q_patch(qp),
					energy(e) {
	}
};

class CCGParticle: public BaseParticle{
protected:
public:
    //strange variables
    number th_b_0=0,beta_b_0=0,xu_bending=0.952319757,xk_bending=1.14813301,kb1=0,kb2=0.80,kt_pref=1;
    CCGParticle();
    virtual ~CCGParticle();
    virtual bool has_bond(BaseParticle *p);
    virtual double return_bfactor(int particleIndex);
    virtual void return_kro(int particleIndex,double *k,double *r0);//return both k and r0 as pointer.
    virtual void add_neighbour(BaseParticle *n, double bfact, double ro); //to add neighbours, affected from base paricle is necessary to call bonded function.
    void set_positions();
    std::vector<int> spring_neighbours;
    std::vector<double>Bfactor,ro;
    double radius,strength=1.0f;
    //Patchy Intereaaction
    std::vector<LR_vector> _base_patches;
    number _sigma;
    std::vector<PatchyBond> bonds;
    // int lockedTo=0; //to have patchy interaction between single pair
    // bool multipatch=true;// if multiple patches could attach or not.
};


#endif //end CCGParticle
// CCGParticle.h
// Subhajit

#ifndef CCGParticle_H
#define CCGParticle_H

#include "BaseParticle.h"
#include <algorithm>

class CCGParticle: public BaseParticle{
protected:
public:
    CCGParticle();
    virtual ~CCGParticle();
    virtual bool has_bond(BaseParticle *p);
    virtual double return_bfactor(int particleIndex);
    virtual void return_kro(int particleIndex,double *k,double *r0);//return both k and r0 as pointer.
    virtual void add_neighbour(BaseParticle *n, double bfact, double ro); //to add neighbours, affected from base paricle is necessary to call bonded function.
    std::vector<int> spring_neighbours;
    std::vector<double>Bfactor,ro;
    double radius,strength=1.0f;
    // int lockedTo=0; //to have patchy interaction between single pair
    // bool multipatch=true;// if multiple patches could attach or not.
};


#endif //end CCGParticle
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
    std::vector<int> spring_neighbours;
    std::vector<double>Bfactor;
};


#endif //end CCGParticle
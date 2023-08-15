//
// Created by jonah on 2/25/21.
//

#ifndef ANM_PARTICLE_H
#define ANM_PARTICLE_H


#include "BaseParticle.h"
//#include <set>

/**
 * @brief A customisable particle. Used by ACInteraction.
 */
class ANMParticle: public BaseParticle {
protected:

public:
    ANMParticle();
    virtual ~ANMParticle();

    virtual bool is_rigid_body() { return false; }

    virtual bool is_bonded(BaseParticle *q);
    virtual void add_bonded_neighbor(BaseParticle *nn);
    std::vector<int> bonded_neighs;
};
#endif //ANM_PARTICLE_H

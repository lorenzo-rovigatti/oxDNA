//
// Created by jonah on 2/25/21.
//

#ifndef GS_PARTICLE_H
#define GS_PARTICLE_H


#include "BaseParticle.h"
//#include <set>

/**
 * @brief A customisable particle. Used by CGDNAInteraction.
 * A literal replica of ANMParticle
 */
class GSParticle: public BaseParticle {
protected:

public:
    GSParticle();
    virtual ~GSParticle();

    virtual bool is_rigid_body() { return false; }

    virtual bool is_bonded(BaseParticle *q);
    virtual void add_bonded_neighbor(BaseParticle *nn);
    std::vector<int> bonded_neighs;
};
#endif //GS_PARTICLE_H

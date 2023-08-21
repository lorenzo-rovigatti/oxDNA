//
// Created by jonah on 2/25/21.
//

#ifndef ANMT_PARTICLE_H
#define ANMT_PARTICLE_H


#include "BaseParticle.h"

/**
 * @brief A customisable particle. Used by ACInteraction.
 */
class ANMTParticle: public BaseParticle {
protected:

public:
    const static LR_vector principal_axis;
    const static LR_vector second_axis;
    const static LR_vector third_axis;

    ANMTParticle();
    virtual ~ANMTParticle();

    virtual bool is_rigid_body() { return true; }

    virtual bool is_bonded(BaseParticle *q);
    virtual void add_bonded_neighbor(BaseParticle *nn);
    std::vector<int> bonded_neighs;
};
#endif //ANMT_PARTICLE_H

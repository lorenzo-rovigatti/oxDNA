//
// Created by jonah,subhajit on feb-2026.
//

#ifndef ANM_PARTICLE_H
#define ANM_PARTICLE_H


#include "BaseParticle.h"

/**
 * @brief A customisable particle. Used by ACInteraction.
 */
class ANMParticle: public BaseParticle {
protected:
    bool _is_rigid_body;

public:
    explicit ANMParticle(bool is_rigid_body = false);
    virtual ~ANMParticle();

    virtual bool is_rigid_body() { return _is_rigid_body; }

    virtual bool is_bonded(BaseParticle *q);
    virtual void add_bonded_neighbor(BaseParticle *nn);
    std::vector<int> bonded_neighs;
    number mass=1.0,invmass=1.0; // Subho
};
#endif //ANM_PARTICLE_H

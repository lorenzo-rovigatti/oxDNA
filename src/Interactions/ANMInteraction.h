//
// Created by jonah on 2/25/21.
//

#ifndef ANM_INTERACTION_H
#define ANM_INTERACTION_H

#include "BaseInteraction.h"
#include "../Particles/ANMParticle.h"

/**
 * @brief Handles (generalised) Anisotropic-Chain interactions between spheres of size .1573 simulation units
 *
 *
 * This interaction is selected with
 * interaction_type = AC
 *
 * Input options:
 *
 * @verbatim
 * none
 */


class ANMInteraction: public BaseInteraction {
protected:

    std::map< std::pair<int, int>, double> _rknot; //eqdist of each bond of psuedobonds
    std::map< std::pair<int, int>, std::pair<char, double> > _potential; //switch to tell lj, FENE or spring as well as strength for each pair of particles
    // switch is currently ignored and just uses spring

    number _rstar, _rc, _b, _sigma;

    inline number _exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    inline number _spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces );
    inline number _repulsive_lj_quart(const LR_vector &r, LR_vector &force, bool update_forces);

public:
    enum {
        SPRING_POTENTIAL = 0,
        EXC_VOL = 1
    };

    ANMInteraction();
    virtual ~ANMInteraction();

    virtual void get_settings(input_file &inp);
    virtual void init();

    virtual void allocate_particles(std::vector<BaseParticle *> &particles);

    virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
    std::map<int, number> masses;
    void load_massfile(std::string &massfile);
    std::string _massfile;

    virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
    virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
    virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);

    virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};



number ANMInteraction::_repulsive_lj_quart(const LR_vector &r, LR_vector &force, bool update_forces) {
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;
    if(rnorm < SQR(_rc)) {
        if(rnorm > SQR(_rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - _rc;
            energy = EXCL_EPS * _b * SQR(SQR(rrc));
            if(update_forces) force = -r * (EXCL_EPS * 4.f * _b * CUB(rrc)/ rmod);
        }
        else {
            number tmp = SQR(_sigma) / rnorm;
            number lj_part = CUB(tmp);
            energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
            if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
        }
    }

    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

    return energy;
}

number ANMInteraction::_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    LR_vector force(0,0,0);
    number energy =  _repulsive_lj_quart(_computed_r, force, update_forces);

    if(update_forces)
    {
        p->force -= force;
        q->force += force;
    }

    return energy;
}

number ANMInteraction::_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {

    std::pair <int,int> keys (std::min(p->index, q->index), std::max(p->index, q->index));

    //Harmonic Spring Potential
    number _k = this->_potential[keys].second; //stiffness of the spring
//    if(p->index == 1 && q->index == 2)
//        printf("eqdist %f spring %f", eqdist, _k);
    number rinsta = _computed_r.module();
    number disp = rinsta - _rknot[keys]; // current distance - eqdistance
    number energy = 0.5 * _k * SQR(disp);

    if (update_forces) {
        LR_vector force(_computed_r);
        force *= (-1.0f * _k ) * disp/rinsta;

        p->force -= force;
        q->force += force;
    }

    return energy;
}

#endif //ANM_INTERACTION_H

//The Below function is currently unused
//    number _repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces);
//
//template<typename number>
//number ANMInteraction::_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces) {
//
//    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
//    number energy = (number) 0;
//    if(rnorm < SQR(_rc)) {
//        if(rnorm > SQR(_rstar)) {
//            number rmod = sqrt(rnorm);
//            number rrc = rmod - _rc;
//            energy = EXCL_EPS * _b * SQR(rrc);
//            if(update_forces) force = -r * (2 * EXCL_EPS * _b * rrc/ rmod);
//        }
//        else {
//            number tmp = SQR(_sigma) / rnorm;
//            number lj_part = tmp * tmp * tmp;
//            energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
//            if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
//        }
//    }
//
//    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;
//
//    return energy;
//}
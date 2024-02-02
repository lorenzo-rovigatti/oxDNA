/*
 * Molecule.cpp
 *
 *  Created on: 21 apr 2020
 *      Author: lorenzo
 */

#include "Molecule.h"

#include "../Boxes/BaseBox.h"
#include "../Forces/ConstantRateTorque.h"
#include "../Forces/ConstantRateForce.h"
#include "../Forces/MovingTrap.h"
#include "../Utilities/ConfigInfo.h"

int Molecule::_current_id = 0;

Molecule::Molecule() : _id(_next_id()) {
}

Molecule::~Molecule() {
}

void Molecule::add_particle(BaseParticle *p) {
    particles.push_back(p);
    _shiftable_dirty = true;
}

void Molecule::normalise() {
    throw oxDNAException("The Molecule::normalise() method has not been implemented yet!");
}

unsigned int Molecule::N() {
    return particles.size();
}

void Molecule::update_com() {
    com = LR_vector(0., 0., 0.);
    for (auto p : particles) {
        com += p->pos;
    }
    com /= N();
}

bool Molecule::shiftable() {
    if (_shiftable_dirty) {
        _is_shiftable = true;
        for (auto p : particles) {
            for (auto force : p->ext_forces) {
                if (typeid(*force) == typeid(MovingTrap) ||
                    typeid(*force) == typeid(ConstantRateTorque) || 
					typeid(*force) == typeid(ConstantRateForce)) {
                    _is_shiftable = false;
                }
            }
        }
        _shiftable_dirty = false;
    }

    return _is_shiftable;
}

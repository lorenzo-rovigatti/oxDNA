//
// Created by jonah on 2/25/21.
//

#include "ANMInteraction.h"

/*
 * ANMInteraction.cpp
 *
 *  Created on: feb-2026
 *      Author: jonah,Subhajit
 */

#include "../Particles/ANMParticle.h"

#include "../Utilities/Utils.h"
#include "ANMTopology.h"

#include <sstream>
#include <unistd.h>
#include <fstream>

// TODO: Get the files loaded with their strings saved in the init



ANMInteraction::ANMInteraction() :
                BaseInteraction() {
    ADD_INTERACTION_TO_MAP(SPRING_POTENTIAL, _spring);
    ADD_INTERACTION_TO_MAP(EXC_VOL, _exc_volume);

    // default mass array
    masses = {{5, 1.f}, {6, 1.f}, {7, 1.f}, {8, 1.f}, {9, 1.f}, {10, 1.f},
            {11, 1.f}, {12, 1.f}, {13, 1.f}, {14, 1.f}, {15, 1.f}, {16, 1.f},
            {17, 1.f}, {18, 1.f}, {19, 1.f}, {20, 1.f}, {21, 1.f}, {22, 1.f},
            {23, 1.f}, {24, 1.f}};
}


ANMInteraction::~ANMInteraction() {

}


void ANMInteraction::get_settings(input_file &inp) {
    BaseInteraction::get_settings(inp);
    char parameterfile[500];
    getInputString(&inp, "parfile", parameterfile, 1);

    //Addition of Reading Parameter File for ANMInteraction Only!
    int key1, key2;
    char potswitch;
    double potential, dist;
    std::string carbons;
    std::fstream parameters;
    parameters.open(parameterfile, std::ios::in);
    getline (parameters,carbons);
    if (parameters.is_open())
    {
        while (parameters >> key1 >> key2 >> dist >> potswitch >> potential)
        {
            std::pair <int, int> lkeys (key1, key2);
            std::pair <char, double> pot (potswitch, potential);
            _rknot[lkeys] = dist;
            _potential[lkeys] = pot;
        }
        parameters.close();
    }
    else
    {
        throw oxDNAException("ParameterFile Could Not Be Opened");
    }

    //Mass File Reading
    if(getInputString(&inp, "massfile", _massfile, 0) == KEY_NOT_FOUND) {
        OX_LOG(Logger::LOG_INFO, "Using Default Mass Array");
    } else {
        load_massfile(_massfile);
    }
}


void ANMInteraction::init() {
    _sqr_rcut = 0.12454; //_rc ^2
    _sigma = 0.35f;
    _rstar = 0.349f;
    _b = 306484596.421f;
    _rc = 0.352894;
}


void ANMInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
    for(auto & particle : particles) {
        particle = new ANMParticle();
    }
}


void ANMInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
    auto topology = ANMTopology::read(this->_topology_filename, "");
    if(topology.N != (int) particles.size()) {
        throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match.");
    }
    allocate_particles(particles);
    ANMTopology::fill_particles(topology, particles, false);
    // as in DNANM, the bonds are the springs of the parameter file (add_bonded_neighbor also fills the affected vectors)
    ANMTopology::add_spring_bonds(_rknot, particles);
    for(auto p : particles) {
        auto *q = static_cast<ANMParticle *>(p);
        q->mass = masses[q->btype];
        q->invmass = 1.f / q->mass;
    }
    *N_strands = topology.N_strands;
}


number ANMInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    number energy;
    if (p->is_bonded(q))
        energy = pair_interaction_bonded(p, q, compute_r, update_forces);
    else
        energy = pair_interaction_nonbonded(p, q, compute_r, update_forces);

    return energy;
}


number ANMInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r && (q != P_VIRTUAL && p != P_VIRTUAL)) {
        _computed_r = q->pos - p->pos;
    }

    number energy = _spring(p,q,compute_r,update_forces);
    energy += _exc_volume(p,q,compute_r,update_forces);
    return energy;
}


number ANMInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r) {
        _computed_r = _box->min_image(p->pos, q->pos);
    }
    if(_computed_r.norm() >= _sqr_rcut) {
        return (number) 0;
    }
    number energy = _exc_volume(p, q, compute_r, update_forces);
    return energy;
}

void ANMInteraction::load_massfile(std::string &filename) {
    std::fstream mass_stream;
    int masstypes;
    mass_stream.open(filename, std::ios::in);
    if(mass_stream.is_open()) {
        int type;
        number mass;
        mass_stream >> masstypes;
        while (mass_stream >> type >> mass) {
            masses[type] = mass;
        }
    } else
        throw oxDNAException("Could Not Load Mass File, Aborting");
}

void ANMInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
}

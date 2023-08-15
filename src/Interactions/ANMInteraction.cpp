//
// Created by jonah on 2/25/21.
//

#include "ANMInteraction.h"

/*
 * ANMInteraction.cpp
 *
 *  Created on: 10/feb/2018
 *      Author: jonah
 */

#include "../Particles/ANMParticle.h"
#include <sstream>
#include <unistd.h>


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
    int N_from_conf = particles.size();
    BaseInteraction::read_topology(N_strands, particles);

    int my_N, my_N_strands;

    char line[5120];
    std::ifstream topology;
    topology.open(this->_topology_filename, std::ios::in);

    if (!topology.good())
        throw oxDNAException("Can't read topology file '%s'. Aborting",
                             this->_topology_filename);

    topology.getline(line, 5120);

    std::stringstream head(line);
    head >> my_N >> my_N_strands;
    if (head.fail()) throw oxDNAException("Problem with header make sure the format is correct for ANM Interaction");

    if(my_N < 0 || my_N < my_N_strands || my_N_strands < 0){
        throw oxDNAException("Problem with header make sure the format is correct for ANM Interaction");
    }

    char aminoacid[256];
    int strand, i = 0;
    while (topology.good()) {
        topology.getline(line, 5120);
        if (strlen(line) == 0 || line[0] == '#')
            continue;
        if (i == N_from_conf)
            throw oxDNAException(
                    "Too many particles found in the topology file (should be %d), aborting",
                    N_from_conf);

        int nside, cside;
        std::stringstream ss(line);
        ss >> strand >> aminoacid >> nside >> cside;
        // This sets the format officially to read from N to C Terminus
        // Same format as most pdb files I've seen

        int x;
        std::set<int> myneighs;
        while(ss.good())
        {
            ss >> x;
            if(x < 0 || x >= N_from_conf)
            {
                throw oxDNAException("Line %d of the topology file has an invalid syntax, neighbor has invalid id",
                                     i + 2);
            }

            myneighs.insert(x);
        }
        auto *p = dynamic_cast<ANMParticle*>(particles[i]);

        if(strlen(aminoacid) == 1) {
            p->type = Utils::decode_aa(aminoacid[0]);
            p->btype = Utils::decode_aa(aminoacid[0]);
            // mass only needs to be set if particle does not have a mass of 1, nucleotides mass =1 by convention
            p->mass = masses[p->btype]; //set mass of each particle
            p->massinverted = 1.f/p->mass; //this needs to be set manually as well
        }

        // add_bonded_neighbor fills affected vector for us
        if (nside < 0)
            p->n3 = P_VIRTUAL;
        else
            p->add_bonded_neighbor(particles[nside]);
        if (cside < 0)
            p->n5 = P_VIRTUAL;
        else
            p->add_bonded_neighbor(particles[cside]);

        for(auto & k : myneighs)
        {
            if(p->index < k) p->add_bonded_neighbor(particles[k]);
        }

        if(p->type == A_INVALID || p->type == P_INVALID)
            throw oxDNAException("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, aminoacid);
        // store the strand id
        // for a design inconsistency, in the topology file
        // strand ids for proteins start from 1, not from 0
        p->strand_id = abs(strand) - 1;
        p->index = i;
        i++;
    }

    if(abs(strand) != my_N_strands) throw oxDNAException("Mismatching strand number in header of Topology file");

    if (i < my_N)
        throw oxDNAException(
                "Not enough particles found in the topology file (should be %d). Aborting",
                my_N);

    topology.close();

    if (my_N != N_from_conf)
        throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");

    *N_strands = my_N_strands;
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

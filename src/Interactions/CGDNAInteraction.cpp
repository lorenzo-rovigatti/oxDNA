/*
 * CGDNAInteraction.cpp
 *
 *  Created on: Aug 13, 2021
 *      Author: jonah
 *  Inherits from the DNA2Interaction Class
 *  Uses DNA2 Model and ANM Protein Model and CG DNA model
 */


#include "CGDNAInteraction.h"
#include <sstream>
#include <fstream>

#include "../Particles/DNANucleotide.h"
#include "../Particles/ANMParticle.h"
#include "../Particles/ANMTParticle.h"

CGDNAInteraction::CGDNAInteraction(bool btp) : DNANMInteraction(btp) { // @suppress("Class members should be properly initialized")
    // TODO: Re-examine These

    //GS-DNA Function Pointers
    ADD_INTERACTION_TO_MAP(GS_EXC_VOL, _gs_exc_volume);
    ADD_INTERACTION_TO_MAP(GS_DNA_EXC_VOL, _gs_dna_exc_volume);
    ADD_INTERACTION_TO_MAP(GS_PRO_EXC_VOL, _gs_pro_exc_volume);

    gs_subtype_num = -1; // gives number of unique gs particles
    _gs_gs_exc_vol_params = nullptr;
    _gs_other_exc_vol_params = nullptr;
    _gs_exc_vol_stiffness = 2.f; // match EXCL_EPS

    topology_order = {-1, -1, -1}; // Number of different particle types

    radii = {{0, 1.f}, {1, 1.f}, {2, 1.f}, {3, 1.f}, {5, 1.f}, {6, 1.f}, //default mass is 1 for everything
              {7, 1.f}, {8, 1.f}, {9, 1.f}, {10, 1.f},{11, 1.f}, {12, 1.f},
              {13, 1.f}, {14, 1.f}, {15, 1.f}, {16, 1.f}, {17, 1.f}, {18, 1.f},
              {19, 1.f}, {20, 1.f}, {21, 1.f}, {22, 1.f}, {23, 1.f}, {24, 1.f}};

}

void CGDNAInteraction::get_settings(input_file &inp){
    this->DNANMInteraction::get_settings(inp);
    // reads par file
    // sets up ANMT parameters if needed
    // reads mass file, the masses array is overwritten in the below code

    if(getInputString(&inp, "massfile", _massfile, 0) == KEY_NOT_FOUND) { // variable mass file
        throw oxDNAException("Mass File Must be provided for CGDNA Interaction");
    }else{
        load_extra_file(_massfile);
    } // CGDNAInteraction must have a mass file defined
}

void CGDNAInteraction::check_input_sanity(std::vector<BaseParticle*> &particles){
    //this->DNA2Interaction::check_input_sanity(particles,N);
    //Need to make own function that checks the input sanity
}

void CGDNAInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
    int n = 0;
    for(auto &t: topology_order) {
        if (t == 0) {
            OX_LOG(Logger::LOG_INFO, "CG Particles Initialized");
            for (int i = n; i < n + ngs; i++) particles[i] = new ANMParticle();
            n += ngs;
        } else if (t == 1) {
            OX_LOG(Logger::LOG_INFO, "Protein Particles Initialized");
            if (_angular) for (int i = n; i < n + npro; i++) particles[i] = new ANMTParticle();
            else for (int i = n; i < n + npro; i++) particles[i] = new ANMParticle();
            n += npro;
        } else if (t == 2) {
            OX_LOG(Logger::LOG_INFO, "DNA Particles Initialized");
            for (int i = n; i < n + ndna; i++) particles[i] = new DNANucleotide(this->_grooving);
            n += ndna;
        }
        // ignore -1 topology_order is initialized with
    }

    assert(n == particles.size());

}

void CGDNAInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
    int my_N, my_N_strands;

    char line[5120];
    std::ifstream topology;
    topology.open(this->_topology_filename, std::ios::in);

    if (!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting",this->_topology_filename);

    topology.getline(line, 5120);
    std::stringstream head(line);

    head >> my_N >> my_N_strands >> ndna >> npro >> ndnas >> npep;
    ngstrands = my_N_strands - ndnas - npep;
    ngs = my_N - ndna - npro;
    if (head.fail()) throw oxDNAException("Problem with header make sure the format is correct for CGDNA Interaction");

    if(my_N_strands < 0 || my_N_strands > my_N || ndna > my_N || ndna < 0 || npro > my_N || npro < 0 || ndnas < 0
            || ndnas > my_N || npep < 0 || npep > my_N || ngstrands < 0 || ngstrands > my_N) {
        throw oxDNAException("Problem with header make sure the format is correct for CGDNA Interaction");
    }

    int strand = 0, i = 0;
    std::string base_tmp;
    std::array<bool, 3> found = {false, false, false}; // specific to models [found_dna, found_protein, found_cg]
    while (topology.good()) {
        topology.getline(line, 5120);
        if (strlen(line) == 0 || line[0] == '#')
            continue;

        std::stringstream ss(line);
        ss >> strand >> base_tmp;
        // CG DNA flag
        // negative strand number and contains digits in the base field "gs0"
        if (strand < 0 && !found[0] && std::any_of(base_tmp.begin(), base_tmp.end(),::isdigit)) {
            topology_order[std::count(found.begin(), found.end(), true)] = 0;
            found[0] = true;
        }
        // Protein Flag
        // negative strand number and no digits in the base field -1 T
        if (strand < 0 && !found[1] && !std::any_of(base_tmp.begin(), base_tmp.end(),::isdigit)) {
            topology_order[std::count(found.begin(), found.end(), true)] = 1;
            found[1] = true;
        }
        // DNA Flag
        if (strand > 0 && !found[2]){ // positive strand number DNA 1 A
            topology_order[std::count(found.begin(), found.end(), true)] = 2;
            found[2] = true;
        }
    }

    topology.seekg(0, std::ifstream::beg); // go back to start of the stream
    topology.getline(line, 5120); // Read the header so we can skip it

    strand = 0, i = 0;
    while (topology.good()) {
        topology.getline(line, 5120);
        if (strlen(line) == 0 || line[0] == '#')
            continue;
        if (i == my_N)
            throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", my_N);

        std::stringstream ss(line);
        ss >> strand;
        if (i == 0) {
             //Must be set prior to allocation of particles
            allocate_particles(particles);
            for (int j = 0; j < my_N; j++) {
                particles[j]->index = j;
                particles[j]->type = P_INVALID;
                particles[j]->strand_id = 0;
            }
        }

        // Amino Acid or... Generic Sphere
        if (strand < 0) {
            BaseParticle *p = particles[i];
            char aminoacid[256];
            int nside, cside;
            int subtype; // used only for GS
            ss >> aminoacid >> nside >> cside;
            int x;
            std::set<int> myneighs;
            std::string tmp(aminoacid);  // conversion to string
            if(std::any_of(tmp.begin(), tmp.end(),::isdigit)){
                //generic sphere topology reading
                tmp.erase(0, 2); // remove first 2 characters -> 'gs', left with subtype integer
                subtype=stoi(tmp);
                if(subtype > gs_subtype_num){
                    gs_subtype_num = subtype;
                }
                p->btype = subtype + 27; // first 27 (0, 26) dna + protein subtypes
                p->type = subtype + 27; // used for interaction matrix calls

                if (nside >= 0) {
                    myneighs.insert(nside);
                    if(_angular) p->n3 = particles[nside];
                }
                if (cside >= 0) {
                    myneighs.insert(cside);
                    if(_angular) p->n5 = particles[cside];
                }

                while (ss.good()) {
                    ss >> x;
                    if (x < 0 || x >= my_N) {
                        throw oxDNAException("Line %d of the topology file has an invalid syntax, neighbor has invalid id", i + 2);
                    }
                    myneighs.insert(x);
                }

                auto *q = dynamic_cast< ANMParticle * > (p);
                for(auto & k : myneighs){
                    if (p->index < k) {
                        q->add_bonded_neighbor(particles[k]);
                    }
                }

                p->mass = masses[p->type];
                p->massinverted = 1.f/p->mass;

                p->strand_id = abs(strand) + ndnas - 1;
                p->index = i;

                i++;

            } else {

                if (nside >= 0) {
                    myneighs.insert(nside);
                    if(_angular) p->n3 = particles[nside];
                }
                if (cside >= 0) {
                    myneighs.insert(cside);
                    if(_angular) p->n5 = particles[cside];
                }
                while (ss.good()) {
                    ss >> x;
                    if (x < 0 || x >= my_N) {
                        throw oxDNAException("Line %d of the topology file has an invalid syntax, neighbor has invalid id", i + 2);
                    }
                    myneighs.insert(x);
                }

                if (_angular) {
                    auto *q = dynamic_cast< ANMTParticle * > (p);
                    for(auto & k : myneighs){
                        if (p->index < k) {
                            q->add_bonded_neighbor(particles[k]);
                        }
                    }

                } else {
                    auto *q = dynamic_cast< ANMParticle * > (p);
                    for(auto & k : myneighs){
                        if (p->index < k) {
                            q->add_bonded_neighbor(particles[k]);
                        }
                    }
                }

                if (strlen(aminoacid) == 1) {
                    p->btype = Utils::decode_aa(aminoacid[0]);
                    p->type = Utils::decode_aa(aminoacid[0]);
                }

                p->mass = masses[p->type];
                p->massinverted = 1.f/p->mass;

                p->strand_id = abs(strand) + ndnas - 1;
                p->index = i;

                i++;
            }
        }
        if (strand > 0) {
            char base[256];
            int tmpn3, tmpn5;
            ss >> base >> tmpn3 >> tmpn5;

            auto *p = dynamic_cast<DNANucleotide *> (particles[i]);

            if (tmpn3 < 0) p->n3 = P_VIRTUAL;
            else p->n3 = particles[tmpn3];
            if (tmpn5 < 0) p->n5 = P_VIRTUAL;
            else p->n5 = particles[tmpn5];


            p->strand_id = strand - 1;

            // the base can be either a char or an integer
            if (strlen(base) == 1) {
                p->btype = Utils::decode_base(base[0]);
                p->type = Utils::decode_base(base[0]);

            } else {
                throw oxDNAException("Only DNA Base Characters Permitted in DNA Strand in Topology");
                //if (atoi(base) > 0) p->type = atoi(base) % 4;
                //else p->type = 3 - ((3 - atoi(base)) % 4);
                //p->btype = atoi(base);
            }

            if (p->type == P_INVALID) throw oxDNAException("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);

            p->mass = 1.f;
            p->massinverted = 1.f;

            p->index = i;
            i++;

            // here we fill the affected vector
            if (p->n3 != P_VIRTUAL) p->affected.emplace_back(ParticlePair(p->n3, p));
            if (p->n5 != P_VIRTUAL) p->affected.emplace_back(ParticlePair(p, p->n5));

        }
        if (strand == 0) throw oxDNAException("No strand 0 should be present please check topology file");

    }

    if (i < my_N)
        throw oxDNAException("Not enough particles found in the topology file (should be %d). Aborting", my_N);

    topology.close();

    if (my_N != (int) particles.size())
        throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");

    *N_strands = my_N_strands;

    number max_rcut = 0;

    // topology must be read for this
    int N = gs_subtype_num;
    int excl_vol_interaction_size = N*N; // triangular matrix
    _gs_gs_exc_vol_params = new number[excl_vol_interaction_size * 4](); // each excl vol needs 5 params per interaction
    for(i = 0; i<excl_vol_interaction_size * 4; i++) _gs_gs_exc_vol_params[i] = 0.f;
    number sigma, rstar, rc, b, interaction_dist;
    int j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i > j)
                continue;
            interaction_dist = radii[27 + i] + radii[27 + j];
            _fill_in_constants(interaction_dist, sigma, rstar, b, rc);
            // for debug
            //printf("interaction dist %.3f \n", interaction_dist);
            //printf("sigma %.3f rstar %.3f b %.3f rc %.3f \n", sigma, rstar, b, rc);
            if(rc > max_rcut) max_rcut = rc;
            // i*gs_subtype_num-(i-1)+j
            _gs_gs_exc_vol_params[4 * (i*N + j)] = sigma;
            _gs_gs_exc_vol_params[4 * (i*N + j) + 1] = rstar;
            _gs_gs_exc_vol_params[4 * (i*N + j) + 2] = b;
            _gs_gs_exc_vol_params[4 * (i*N + j) + 3] = rc;

            _gs_gs_exc_vol_params[4 * (j*N + i)] = sigma;
            _gs_gs_exc_vol_params[4 * (j*N + i) + 1] = rstar;
            _gs_gs_exc_vol_params[4 * (j*N + i) + 2] = b;
            _gs_gs_exc_vol_params[4 * (j*N + i) + 3] = rc;
        }
    }

    int excl_vol_pro_dna_size = N * 3;
    _gs_other_exc_vol_params = new number[excl_vol_pro_dna_size * 4]();
    for(i = 0; i<excl_vol_pro_dna_size * 4; i++) _gs_other_exc_vol_params[i] = 0.f;
    // fill in protein and dna excl volume
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            if (j == 0) {
                _fill_in_constants(radii[i+27] + _pro_sigma / 2, sigma, rstar, b, rc); // protein 0
            } else if (j == 1) {
                _fill_in_constants(radii[i+27] + _pro_base_sigma - .175f, sigma, rstar, b, rc); // base 1
            } else if (j == 2) {
                _fill_in_constants(radii[i+27] + _pro_backbone_sigma - .175f, sigma, rstar, b, rc); // backbone 2
            } else {
                throw oxDNAException("Problem filling in Excluded Volume Constants");
            }
            if(rc > max_rcut) max_rcut = rc;
            // for debug
            //printf("sigma %.3f rstar %.3f b %.3f rc %.3f \n", sigma, rstar, b, rc);
            _gs_other_exc_vol_params[4 * (3 * i + j)] = sigma;
            _gs_other_exc_vol_params[4 * (3 * i + j) + 1] = rstar;
            _gs_other_exc_vol_params[4 * (3 * i + j) + 2] = b;
            _gs_other_exc_vol_params[4 * (3 * i + j) + 3] = rc;
        }
    }

    // max interaction distance
    _rcut = max_rcut;
}

// from type (0-4) DNA, (5-26) Protein, (27+) CGDNA
// return index 0 for dna, 1 for protein and 2 for cgdna
int CGDNAInteraction::particle_id(int type){
    return (type <= 4) ? 0: (type >= 27) ? 2 : 1;
}

number CGDNAInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    int pid = particle_id(p->type);
    int qid = particle_id(q->type);
    if(pid == qid){
        if(p->is_bonded(q)) return pair_interaction_bonded(p, q, compute_r, update_forces);
        else return pair_interaction_nonbonded(p, q, compute_r, update_forces);
    } else {
        return pair_interaction_nonbonded(p, q, compute_r, update_forces);
    }
}

number CGDNAInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r)
        if (q != P_VIRTUAL && p != P_VIRTUAL)
            _computed_r = this->_box->min_image(p->pos, q->pos);

    int pid = particle_id(p->type);
    if(p->is_bonded(q)) return (this->*_interaction_matrix_bonded[pid])(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_gs_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    // get params first
    number energy = _protein_spring(p, q, compute_r, update_forces);
    energy += _gs_exc_volume(p, q, compute_r, update_forces);
    return energy;
}

number CGDNAInteraction::_pro_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    if(_angular){
        if (!p->is_bonded(q)) return 0.f;
        number energy = _protein_spring(p, q, compute_r, update_forces);
        energy += _protein_exc_volume(p, q, compute_r, update_forces);
        if (q->index - p->index == 1) energy += _protein_ang_pot(p, q, compute_r, update_forces);
        return energy;
    } else {
        if (!p->is_bonded(q)) return 0.f;
        number energy = _protein_spring(p, q, compute_r, update_forces);
        energy += _protein_exc_volume(p, q, compute_r, update_forces);
        return energy;
    }
}

number CGDNAInteraction::_dna_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    if(!this->_check_bonded_neighbour(&p, &q, compute_r)) return (number) 0;
    number energy = _backbone(p,q,compute_r,update_forces);
    energy += _bonded_excluded_volume(p,q,compute_r,update_forces);
    energy += _stacking(p,q,compute_r,update_forces);
    return energy;
}

number CGDNAInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if (compute_r)
        _computed_r = this->_box->min_image(p->pos, q->pos);

    int pid = particle_id(p->type);
    int qid = particle_id(q->type);
    return (this->*_interaction_matrix_nonbonded[3*pid + qid])(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_dna_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    number rnorm = _computed_r.norm();
    if (rnorm >= _sqr_rcut) return (number) 0.f;
    number energy = _nonbonded_excluded_volume(p, q, compute_r, update_forces);
    energy += _hydrogen_bonding(p, q, compute_r, update_forces);
    energy += _cross_stacking(p, q, compute_r, update_forces);
    energy += _coaxial_stacking(p, q, compute_r, update_forces);
    energy += _debye_huckel(p, q, compute_r, update_forces);
    return energy;
}

number CGDNAInteraction::_pro_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    number rnorm = _computed_r.norm();
    if (rnorm >= _pro_sqr_rcut) return (number) 0.f;
    else return _protein_exc_volume(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_gs_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    return _gs_exc_volume(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_dna_pro_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    number rnorm = _computed_r.norm();
    if (rnorm >= _pro_dna_sqr_rcut) return (number) 0.f;
    return _protein_dna_exc_volume(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_dna_gs_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    return _gs_dna_exc_volume(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_pro_gs_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    return _gs_pro_exc_volume(p, q, compute_r, update_forces);
}

number CGDNAInteraction::_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness) {
    // this is a bit faster than calling r.norm()
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;

    if(rnorm < SQR(rcut)) {
        if(rnorm > SQR(rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - rcut;
            energy = stiffness * b * SQR(rrc);
            if(update_forces) force = -r * (4 * stiffness * b * CUB(rrc) / rmod);
        }
        else {
            number tmp = SQR(sigma) / rnorm;
            number lj_part = tmp * tmp * tmp;
            energy = 4 * stiffness * (SQR(lj_part) - lj_part);
            if(update_forces) force = -r * (24 * stiffness * (lj_part - 2*SQR(lj_part)) / rnorm);
        }
    }

    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

    return energy;
}

number CGDNAInteraction::_repulsive_lj_quart(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness) {
    // this is a bit faster than calling r.norm()
    //changed to a quartic form
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;
    if(rnorm < SQR(rcut)) {
        if(rnorm > SQR(rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - rcut;
            energy = EXCL_EPS * b * SQR(SQR(rrc));
            if(update_forces) force = -r * (4 * EXCL_EPS * b * CUB(rrc)/ rmod);
        }
        else {
            number tmp = SQR(sigma) / rnorm;
            number lj_part = tmp * tmp * tmp;
            energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
            if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
        }
    }

    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

    return energy;
}

number CGDNAInteraction::_gs_dna_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    BaseParticle *gs;
    BaseParticle *nuc;

    LR_vector force(0, 0, 0);
    LR_vector rcenter = _computed_r;

    int pid = particle_id(p->type);

    if(pid == 0)
    {
        gs = q;
        nuc = p;
    } else {
        rcenter = -rcenter;
        gs = p;
        nuc = q;
    }

    LR_vector r_to_back = rcenter  - nuc->int_centers[DNANucleotide::BACK];
    LR_vector r_to_base = rcenter  - nuc->int_centers[DNANucleotide::BASE];

    LR_vector torquenuc(0,0,0);
    auto energy = (number) 0.f;

    int i = gs->type - 27;

    number sigma, rstar, b, rc;

    // backbone parameters
    sigma = _gs_other_exc_vol_params[4 * (3*i + 2)];
    rstar = _gs_other_exc_vol_params[4 * (3*i + 2) + 1];
    b = _gs_other_exc_vol_params[4 * (3*i + 2) + 2];
    rc = _gs_other_exc_vol_params[4 * (3*i + 2) + 3];

    energy += _repulsive_lj(r_to_back, force, update_forces, sigma, b, rstar,rc,_gs_exc_vol_stiffness);
    if (update_forces) {
        torquenuc = nuc->int_centers[DNANucleotide::BACK].cross(-force);
        nuc->torque += nuc->orientationT * torquenuc;
        nuc->force -= force;
        gs->force += force;
    }

    // base parameters
    sigma = _gs_other_exc_vol_params[4 * (3 * i + 1)];
    rstar = _gs_other_exc_vol_params[4 * (3 * i + 1) + 1];
    b = _gs_other_exc_vol_params[4 * (3 * i + 1) + 2];
    rc = _gs_other_exc_vol_params[4 * (3 * i + 1) + 3];

    energy += _repulsive_lj(r_to_base, force, update_forces, sigma, b, rstar, rc, _gs_exc_vol_stiffness);
    if(update_forces) {
        torquenuc = nuc->int_centers[DNANucleotide::BASE].cross(-force);
        nuc->torque += nuc->orientationT * torquenuc;
        nuc->force -= force;
            gs->force += force;
    }

    return energy;
}

number CGDNAInteraction::_gs_pro_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    int i;
    p->type < q->type ? i=q->type: i=p->type; // highest type number is gs, assign type to i
    i -= 27;

    number sigma = _gs_other_exc_vol_params[4 * (3*i + 0)];
    number rstar = _gs_other_exc_vol_params[4 * (3*i + 0) + 1];
    number b = _gs_other_exc_vol_params[4 * (3*i + 0) + 2];
    number rc = _gs_other_exc_vol_params[4 * (3*i + 0) + 3];

    LR_vector force(0, 0, 0);
    number energy = _repulsive_lj(_computed_r, force, update_forces, sigma, b, rstar, rc, _gs_exc_vol_stiffness);
    if(update_forces)
    {
        p->force -= force;
        q->force += force;
    }
    return energy;
}

number CGDNAInteraction::_gs_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    int i = p->type-27; // ignores the first 27 types defined by dna bases and amino acids
    int j = q->type-27;
    int &N = gs_subtype_num;

    number sigma = _gs_gs_exc_vol_params[4 * (i*N+ j)];
    number rstar = _gs_gs_exc_vol_params[4 * (i*N+ j)+1];
    number b = _gs_gs_exc_vol_params[4 * (i*N+ j)+2];
    number rc = _gs_gs_exc_vol_params[4 * (i*N+ j)+3];

    LR_vector force(0, 0, 0);
    number energy = _repulsive_lj(_computed_r, force, update_forces, sigma, b, rstar, rc, _gs_exc_vol_stiffness);
    if(update_forces)
    {
        p->force -= force;
        q->force += force;
    }
    return energy;
}

void CGDNAInteraction::init() {
    DNANMInteraction::init();
    ngs = 0, ngstrands = 0, npep = 0;
    // N(N+1)/2 possible pair interactions for N subtypes of gs particles
    //int N = gs_subtype_num;
}

void CGDNAInteraction::load_extra_file(std::string &filename) {
    std::fstream mass_stream;
    int masstypes;
    mass_stream.open(filename, std::ios::in);
    if(mass_stream.is_open()) {
        int type;
        number mass;
        number radius;
        mass_stream >> masstypes;
        masses.clear(); // remove default masses
        radii.clear(); // remove default radii
        while (mass_stream >> type >> mass >> radius) {
            masses[type] = (number) mass;
            radii[type] = (number) radius;
        }
    } else
        throw oxDNAException("Could Not Load Mass File, Aborting");
}

CGDNAInteraction::~CGDNAInteraction() {
    // clean up our pointers with allocated memory
    delete[] _gs_gs_exc_vol_params;
    delete[] _gs_other_exc_vol_params;
}






/*
 * DNANMInteraction.cpp
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 *  Inherits from the DNA2Interaction Class
 *  Uses DNA2 Model and ANM Protein Model
 */


#include "DNANMInteraction.h"
#include <sstream>
#include <fstream>

#include "../Particles/DNANucleotide.h"
#include "../Particles/ANMParticle.h"
#include "../Particles/ANMTParticle.h"


DNANMInteraction::DNANMInteraction(bool btp) : DNA2Interaction() {

    _angular = btp;
    //Protein Methods Function Pointers
    ADD_INTERACTION_TO_MAP(SPRING, _protein_spring);
    ADD_INTERACTION_TO_MAP(PRO_EXC_VOL, _protein_exc_volume);
    if(_angular)
        ADD_INTERACTION_TO_MAP(PRO_ANG_POT, _protein_ang_pot);

    //Protein-DNA Function Pointers
    ADD_INTERACTION_TO_MAP(PRO_DNA_EXC_VOL, _protein_dna_exc_volume);


    //DNA Methods Function Pointers
    ADD_INTERACTION_TO_MAP(BACKBONE, _backbone);
    ADD_INTERACTION_TO_MAP(COAXIAL_STACKING, _coaxial_stacking);
    ADD_INTERACTION_TO_MAP(CROSS_STACKING, _cross_stacking);
    ADD_INTERACTION_TO_MAP(BONDED_EXCLUDED_VOLUME, _bonded_excluded_volume);
    ADD_INTERACTION_TO_MAP(STACKING, _stacking);
    ADD_INTERACTION_TO_MAP(HYDROGEN_BONDING, _hydrogen_bonding);
    ADD_INTERACTION_TO_MAP(NONBONDED_EXCLUDED_VOLUME, _nonbonded_excluded_volume);
    ADD_INTERACTION_TO_MAP(DEBYE_HUCKEL, _debye_huckel);

    _parameter_kbkt = false;

    masses = {{0, 1.f}, {1, 1.f}, {2, 1.f}, {3, 1.f}, {5, 1.f}, {6, 1.f}, //default mass is 1 for everything
              {7, 1.f}, {8, 1.f}, {9, 1.f}, {10, 1.f},{11, 1.f}, {12, 1.f},
              {13, 1.f}, {14, 1.f}, {15, 1.f}, {16, 1.f}, {17, 1.f}, {18, 1.f},
              {19, 1.f}, {20, 1.f}, {21, 1.f}, {22, 1.f}, {23, 1.f}, {24, 1.f}};
}


void DNANMInteraction::get_settings(input_file &inp){
    this->DNA2Interaction::get_settings(inp);

    if(_angular){
        float kb_tmp = -1.f;
        float kt_tmp = -1.f;
        getInputString(&inp, "parfile", _parameterfile, 1);
        getInputFloat(&inp, "bending_k", &kb_tmp, 0);
        _k_bend = (number) kb_tmp;

        getInputFloat(&inp, "torsion_k", &kt_tmp, 0);
        _k_tor = (number) kt_tmp;

        if((_k_bend < 0 && kt_tmp >= 0) || (_k_tor < 0 && _k_bend >= 0)){
            throw oxDNAException("Bending & Torsion Constants Must BOTH be set globally or in parameter file");
        } else if(_k_bend >= 0 && _k_tor >= 0){
            OX_LOG(Logger::LOG_INFO, "Using Global kb/kt values set in Input File");
        } else if(_k_bend < 0 && _k_tor < 0){
            OX_LOG(Logger::LOG_INFO, "No Global kb/kt values set in Input File, Attempting to read from Parameter File");
            _parameter_kbkt = true;
        }
    }

    getInputString(&inp, "parfile", _parameterfile, 0); //parameter file

    if(getInputString(&inp, "massfile", _massfile, 0) == KEY_NOT_FOUND) { // variable mass file
        OX_LOG(Logger::LOG_INFO, "Using Default Masses"); // declared in constructor
    } else {
        OX_LOG(Logger::LOG_INFO, "Using Provided Massfile");
        load_massfile(_massfile);
    }
}


void DNANMInteraction::read_parameter_file(std::vector<BaseParticle*> &particles){
    //Addition of Reading Parameter File
    auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
        if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
        if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
        if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
        if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
        if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
    };

    //Checkers as Lambdas
    auto valid_angles = [](double a, double b, double c, double d)
    {
        double anglemin = std::min({a, b, c, d});
        double anglemax = std::max({a, b, c, d});
        if (anglemin < -1.0 || anglemax > 1.0){
            throw oxDNAException("Cos of Angle in Parameter File not in Valid bounds");
        }
    };

    char n[5] = "none";
    if(strcmp(_parameterfile, n) != 0 && strcmp(_parameterfile, "") != 0) {
        int key1, key2;
        char potswitch;
        double potential, dist;
        double a0, b0, c0, d0;
        int N;
        std::fstream parameters;
        parameters.open(_parameterfile, std::ios::in);
        parameters >> N;
        int spring_connection_num = 0;
        if (parameters.is_open())
        {
            while (parameters >> key1 >> key2 >> dist >> potswitch >> potential)
            {
                valid_spring_params(N, key1, key2, dist, potswitch, potential);
                if (_angular) {
                    auto *q = dynamic_cast< ANMTParticle * > (particles[key1]);
                    if (key1 < key2) {
                        q->add_bonded_neighbor(particles[key2]);
                    }
                } else {
                    auto *q = dynamic_cast< ANMParticle * > (particles[key1]);
                    if (key1 < key2) {
                        q->add_bonded_neighbor(particles[key2]);
                    }
                }

                spring_connection_num += 1;

                //If DNACT
                if(_angular) {
                    if (key2 - key1 == 1) {
                        try {
                            parameters >> a0 >> b0 >> c0 >> d0;
                            valid_angles(a0, b0, c0, d0);
                            if (_parameter_kbkt) {
                                parameters >> _k_bend >> _k_tor;
                                if (_k_bend < 0 || _k_tor < 0)
                                    throw oxDNAException("Invalid pairwise kb/kt Value Declared in Parameter File");
                                std::pair<double, double> ang_constants(_k_bend, _k_tor);
                                _ang_stiff[key1] = ang_constants;
                            }
                            std::vector<double> angles{a0, b0, c0, d0};
                            _ang_vals[key1] = angles;
                        } catch(...){
                            continue;
                        }
                    }
                }
                std::pair <int, int> lkeys (key1, key2);
                std::pair <char, double> pot (potswitch, potential);
                _rknot[lkeys] = dist;
                _potential[lkeys] = pot;
            }
        }
        else
        {
            throw oxDNAException("ParameterFile Could Not Be Opened");
        }
        parameters.close();
        if (_angular) {
            if (spring_connection_num == 1 && N > 2 && _parameter_kbkt)
                throw oxDNAException(
                        "Invalid Parameter File Format, pairwise kb and kt values incorrectly formatted or missing");
            if (spring_connection_num == 1 && N > 2 && !_parameter_kbkt)
                throw oxDNAException(
                        "Invalid Parameter File Format, global kb and kt values have been set in Inputfile");
        } else if (spring_connection_num == 1 && N > 2) throw oxDNAException("Invalid Parameter File Format, cannot use a DNACT Parameter File");
    } else {
        OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");
    }
}


void DNANMInteraction::check_input_sanity(std::vector<BaseParticle*> &particles){
    //this->DNA2Interaction::check_input_sanity(particles,N);
    //Need to make own function that checks the input sanity
}


void DNANMInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
    if (ndna==0 || ndnas==0) {
        OX_LOG(Logger::LOG_INFO, "No DNA Particles Specified, Continuing with just Protein Particles");
        if (_angular) for (int i = 0; i < npro; i++) particles[i] = new ANMTParticle();
        else for (int i = 0; i < npro; i++) particles[i] = new ANMParticle();

    } else if (npro == 0) {
        OX_LOG(Logger::LOG_INFO, "No Protein Particles Specified, Continuing with just DNA Particles");
        for (int i = 0; i < ndna; i++) particles[i] = new DNANucleotide(this->_grooving);

    } else {
        if (_firststrand > 0){
            for (int i = 0; i < ndna; i++) particles[i] = new DNANucleotide(this->_grooving);
            // Protein
            if (_angular) for (uint i = ndna; i < particles.size(); i++) particles[i] = new ANMTParticle();
            else for (uint i = ndna; i < particles.size(); i++) particles[i] = new ANMParticle();
        } else {
            // Protein
            if (_angular) for (int i = 0; i < npro; i++) particles[i] = new ANMTParticle();
            else for (int i = 0; i < npro; i++) particles[i] = new ANMParticle();

            for (uint i = npro; i < particles.size(); i++) particles[i] = new DNANucleotide(this->_grooving);
        }
    }
}


void DNANMInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
    int my_N, my_N_strands;

    char line[5120];
    std::ifstream topology;
    topology.open(this->_topology_filename, std::ios::in);

    if (!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting",this->_topology_filename);

    topology.getline(line, 5120);
    std::stringstream head(line);

    head >> my_N >> my_N_strands >> ndna >> npro >>ndnas;
    if (head.fail()) throw oxDNAException("Problem with header make sure the format is correct for DNANM Interaction");

    if(my_N_strands < 0 || my_N_strands > my_N || ndna > my_N || ndna < 0 || npro > my_N || npro < 0 || ndnas < 0 || ndnas > my_N) {
        throw oxDNAException("Problem with header make sure the format is correct for DNANM Interaction");
    }


    int strand, i = 0;
    while (topology.good()) {
        topology.getline(line, 5120);
        if (strlen(line) == 0 || line[0] == '#')
            continue;
        if (i == my_N)
            throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", my_N);

        std::stringstream ss(line);
        ss >> strand;
        if (i == 0) {
            _firststrand = strand; //Must be set prior to allocation of particles
            allocate_particles(particles);
            for (int j = 0; j < my_N; j++) {
                particles[j]->index = j;
                particles[j]->type = P_INVALID;
                particles[j]->strand_id = 0;
            }
        }

        // Amino Acid
        if (strand < 0) {
            BaseParticle *p = particles[i];
            char aminoacid[256];
            int nside, cside;
            ss >> aminoacid >> nside >> cside;

            int x;
            std::set<int> myneighs;
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

            if (strlen(aminoacid) == 1) {
                p->type = Utils::decode_aa(aminoacid[0]);
                p->btype = Utils::decode_aa(aminoacid[0]);; // btype of 1 means AA
            }

            p->mass = masses[p->type];
            p->massinverted = 1.f/p->mass;

            p->strand_id = abs(strand) + ndnas - 1;
            p->index = i;

            i++;
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
                p->type = Utils::decode_base(base[0]);
                p->btype = Utils::decode_base(base[0]);

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

    // read parameter file
    read_parameter_file(particles);
}


number DNANMInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    int interaction_type = get_id(p->btype) + get_id(q->btype);

    if (interaction_type == 0 || interaction_type == 2){
        if(p->is_bonded(q)) return pair_interaction_bonded(p, q, compute_r, update_forces);
        else return pair_interaction_nonbonded(p, q, compute_r, update_forces);
    }
    if (interaction_type == 1) return this->pair_interaction_nonbonded(p, q, compute_r, update_forces);

    return 0.f;
}


number DNANMInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r)
        if (q != P_VIRTUAL && p != P_VIRTUAL)
            _computed_r = this->_box->min_image(p->pos, q->pos);

    int interaction_type = get_id(p->btype) + get_id(q->btype);
    if (interaction_type == 0){ // dna-dna
        if(!this->_check_bonded_neighbour(&p, &q, compute_r)) return (number) 0;
        number energy = _backbone(p,q,compute_r,update_forces);
        energy += _bonded_excluded_volume(p,q,compute_r,update_forces);
        energy += _stacking(p,q,compute_r,update_forces);
        return energy;
    } else if (interaction_type == 2){ // protein-protein
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
    } else
        return 0.f;

    // Expect Topology to not contain Bonds b/t protein & DNA
}


number DNANMInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if (compute_r)
        _computed_r = this->_box->min_image(p->pos, q->pos);

    number rnorm = _computed_r.norm();

    int interaction_type = get_id(p->btype) + get_id(q->btype);

    if (interaction_type == 0) { //DNA-DNA Interaction
        if (rnorm >= _sqr_rcut) return (number) 0.f;
        number energy = _nonbonded_excluded_volume(p, q, compute_r, update_forces);
        energy += _hydrogen_bonding(p, q, compute_r, update_forces);
        energy += _cross_stacking(p, q, compute_r, update_forces);
        energy += _coaxial_stacking(p, q, compute_r, update_forces);
        energy += _debye_huckel(p, q, compute_r, update_forces);
        return energy;
    } else if (interaction_type == 2) { // protein-protein
        if (rnorm >= _pro_sqr_rcut) return (number) 0.f;
        number energy = _protein_exc_volume(p, q, compute_r, update_forces);
        return energy;
//        return 0.f;
    }else if(interaction_type == 1) { //protein-dna
        if (rnorm >= _pro_dna_sqr_rcut) return (number) 0.f;
        number energy = _protein_dna_exc_volume(p, q, compute_r, update_forces);
        return energy;
    }

    return 0.f;
}


number DNANMInteraction::_protein_dna_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces)
{
    BaseParticle *protein;
    BaseParticle *nuc;

    LR_vector force(0, 0, 0);
    LR_vector rcenter = _computed_r;
    int pid = get_id(p->btype);
    int qid = get_id(q->btype);

    if(pid == 0 && qid == 1)
    {
        //rcenter = -rcenter;
        protein = q;
        nuc = p;
    }
    else if (pid == 1 && qid == 0)
    {
        rcenter = -rcenter;
        protein = p;
        nuc = q;
    }
    else
        return 0.f;

    LR_vector r_to_back = rcenter  - nuc->int_centers[DNANucleotide::BACK];
    LR_vector r_to_base = rcenter  - nuc->int_centers[DNANucleotide::BASE];

    LR_vector torquenuc(0,0,0);
    auto energy = (number) 0.f;


    if(r_to_back.norm() < _pro_backbone_sqr_rcut) {
        energy = _protein_dna_repulsive_lj(r_to_back, force, update_forces, _pro_backbone_sigma, _pro_backbone_b, _pro_backbone_rstar,_pro_backbone_rcut,_pro_backbone_stiffness);
        //printf("back-pro %d %d %f\n",p->index,q->index,energy);
        if (update_forces) {
            torquenuc = nuc->int_centers[DNANucleotide::BACK].cross(-force);
            nuc->torque += nuc->orientationT * torquenuc;
            nuc->force -= force;
            protein->force += force;
            if(_angular){
                r_to_back.normalize();
                protein->torque += p->orientationT*(-r_to_back*_pro_sigma*0.5f).cross(force); // dr Point of contact on protein particle relative to COM of protein particle
            }
        }
    }

    if(r_to_base.norm() < _pro_base_sqr_rcut){
//        printf("pro %d dna %d\n", protein->index, nuc->index);
        energy += _protein_dna_repulsive_lj(r_to_base, force, update_forces, _pro_base_sigma, _pro_base_b, _pro_base_rstar, _pro_base_rcut, _pro_base_stiffness);
        if(update_forces) {
            torquenuc = nuc->int_centers[DNANucleotide::BASE].cross(-force);
            nuc->torque += nuc->orientationT * torquenuc;
            nuc->force -= force;
            protein->force += force;
            if(_angular){
                r_to_base.normalize();
                protein->torque += p->orientationT*(-r_to_base*_pro_sigma*0.5f).cross(force); // dr Point of contact on protein particle relative to COM of protein particle
            }
        }
    }

    return energy;
}


number DNANMInteraction::_protein_dna_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness) {
    // this is a bit faster than calling r.norm()
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;

    if(rnorm < SQR(rcut)) {
        if(rnorm > SQR(rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - rcut;
            energy = stiffness * b * SQR(SQR(rrc));
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


void DNANMInteraction::init() {
    DNA2Interaction::init();
    ndna=0, npro=0, ndnas =0;
    //let's try this
    _pro_backbone_sigma = 0.57f;
    _pro_backbone_rstar= 0.569f;
    _pro_backbone_b = 178699253.5f;
    _pro_backbone_rcut = 0.572934f;
    _pro_backbone_stiffness = 1.0f;
    _pro_backbone_sqr_rcut = 0.3283f;
    //Base-Protein Excluded Volume Parameters
    _pro_base_sigma = 0.36f;
    _pro_base_rstar= 0.359f;
    _pro_base_b = 296866090.f;
    _pro_base_rcut = 0.362897f;
    _pro_base_stiffness = 1.0f;
    _pro_base_sqr_rcut = 0.1317f;
    //Protein-Protein Excluded Volume Parameters
    _pro_sigma = 0.35f;
    _pro_rstar = 0.349f;
    _pro_b = 306484596.421f;
    _pro_rcut = 0.352894;
    _pro_sqr_rcut = 0.12454f; //_rc ^2

    _pro_dna_sqr_rcut = 3.0625f; //placeholder test value 1.75^2
}

//Functions almost Identical to those in ANMInteraction
number DNANMInteraction::_protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces) {
    // this is a bit faster than calling r.norm()
    //changed to a quartic form
    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
    number energy = (number) 0;
    if(rnorm < SQR(_pro_rcut)) {
        if(rnorm > SQR(_pro_rstar)) {
            number rmod = sqrt(rnorm);
            number rrc = rmod - _pro_rcut;
            energy = EXCL_EPS * _pro_b * SQR(SQR(rrc));
            if(update_forces) force = -r * (4 * EXCL_EPS * _pro_b * CUB(rrc)/ rmod);
        }
        else {
            number tmp = SQR(_pro_sigma) / rnorm;
            number lj_part = tmp * tmp * tmp;
            energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
            if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
        }
    }

    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

    return energy;
}


number DNANMInteraction::_protein_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    LR_vector force(0,0,0);
    number energy =  _protein_repulsive_lj(_computed_r, force, update_forces);

    if(update_forces)
    {
        p->force -= force;
        q->force += force;
        if(_angular){
            p->torque += p->orientationT*(_computed_r*0.5).cross(-force);
            q->torque += q->orientationT*(-_computed_r*0.5).cross(force);
        }
    }

    return energy;
}


number DNANMInteraction::_protein_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {

    std::pair <int,int> keys (std::min(p->index, q->index), std::max(p->index, q->index));

    //Harmonic Spring Potential
    number _k = _potential[keys].second; //stiffness of the spring

    number rinsta = _computed_r.module(); // distance b/t p and q
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


number DNANMInteraction::_protein_ang_pot(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    // Get Angular Parameters
    std::vector<double> &ang_params = _ang_vals[p->index];
    double &a0 = ang_params[0];
    double &b0 = ang_params[1];
    double &c0 = ang_params[2];
    double &d0 = ang_params[3];

    if (_parameter_kbkt) { //Uses array with per particle kb and kt
        _k_bend = this->_ang_stiff[p->index].first;
        _k_tor = this->_ang_stiff[p->index].second;
    } // If False, the global kb and kt will be used


    LR_vector rij_unit = _computed_r;
    rij_unit.normalize();
    LR_vector rji_unit = rij_unit * -1.f;

    LR_vector &a1 = p->orientationT.v1;
    LR_vector &b1 = q->orientationT.v1;
    LR_vector &a3 = p->orientationT.v3;
    LR_vector &b3 = q->orientationT.v3;


    double o1 = rij_unit * a1 - a0;
    double o2 = rji_unit * b1 - b0;
    double o3 = a1 * b1 - c0;
    double o4 = a3 * b3 - d0;

    //Torsion and Bending
    number energy = _k_bend * 0.5f * (SQR(o1) + SQR(o2)) + _k_tor * 0.5f * (SQR(o3) + SQR(o4));

    if (update_forces) {

        LR_vector force = -(rji_unit.cross(rij_unit.cross(a1)) * _k_bend * o1 -
                                    rji_unit.cross(rij_unit.cross(b1)) * _k_bend * o2) / _computed_r.module();

        p->force -= force;
        q->force += force;

        LR_vector ta = rij_unit.cross(a1) * o1 * _k_bend;
        LR_vector tb = rji_unit.cross(b1) * o2 * _k_bend;

        p->torque += p->orientationT * ta;
        q->torque += q->orientationT * tb;

        LR_vector torsion = a1.cross(b1) * o3 * _k_tor + a3.cross(b3) * o4 * _k_tor;

        p->torque += p->orientationT * -torsion;
        q->torque += q->orientationT * torsion;

        //For Debugging, very helpful for CUDA comparisons
        /*LR_vector<number> TA = p->orientationT * (ta - torsion);
        LR_vector<number> TB = q->orientationT * (tb + torsion);

        if (p->index ==104){
            printf("p2 Angular F %.7f %.7f %.7f TA %.7f %.7f %.7f\n", -force.x, -force.y, -force.z, TA.x, TA.y, TA.z);
        } else if (q->index ==104){
            printf("q2 Angular F %.7f %.7f %.7f TB %.7f %.7f %.7f\n", force.x, force.y, force.z, TB.x, TB.y, TB.z);
        }*/

    }

    return energy;
}


void DNANMInteraction::load_massfile(std::string &filename) {
    std::fstream mass_stream;
    int masstypes;
    mass_stream.open(filename, std::ios::in);
    if(mass_stream.is_open()) {
        int type;
        number mass;
        mass_stream >> masstypes;
        masses.clear(); // remove default masses
        while (mass_stream >> type >> mass) {
            masses[type] = (number) mass;
        }
    } else
        throw oxDNAException("Could Not Load Mass File, Aborting");
}

int DNANMInteraction::get_id(int btype){
    // takes btype return whether dna or protein
    // 0 is dna
    // >=5 is protein
    return (btype <= 4) ? 0: 1;
};

DNANMInteraction::~DNANMInteraction() = default;






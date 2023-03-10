
/**
 * DNANMwithRNANMInteraction.cpp
 *
 *  Created on: Jan 16, 2023
 *      Author: jonah (+eryk)
 * 
 * (just the DNANM interaction without templates)
 *
 * 
 * 
 **/

//(for some basic debugging OX_LOG(Logger::LOG_INFO, "We made it this far.");)
#include "DNANMwithRNANMInteraction.h"
#include <sstream>
#include <fstream>
#include <utility>
#include <cfloat>

#include "../Particles/DNANucleotide.h"
#include "../Particles/RNANucleotide.h"
#include "../Particles/ACParticle.h"
#include "rna_model.h"


DNANMwithRNANMInteraction::DNANMwithRNANMInteraction() : DNA2withRNA2Interaction() { // @suppress("Class members should be properly initialized")


    //these cause errors, fix later
    //ADD_INTERACTION_TO_MAP(SPRING, _protein_spring);
    //ADD_INTERACTION_TO_MAP(PRO_EXC_VOL, _protein_exc_volume);
    //ADD_INTERACTION_TO_MAP(PRO_DNA_EXC_VOL, _protein_dna_exc_volume);
}

DNANMwithRNANMInteraction::~DNANMwithRNANMInteraction() {
}

void DNANMwithRNANMInteraction::get_settings(input_file &inp){

	DNA2withRNA2Interaction::get_settings(inp);
	getInputString(&inp, "parfile", _parameterfile, 0);
	//Addition of Reading Parameter File
    char n[5] = "none";

    auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
        if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
        if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
        if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
        if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
        if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
    };

    if(strcmp(_parameterfile, n) != 0) {
        int key1, key2;
        char potswitch;
        double potential, dist;
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
                spring_connection_num += 1;
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
        if(spring_connection_num == 1 && N > 2) throw oxDNAException("Invalid Parameter File Format, cannot use a DNACT Parameter File");
    } else {
        OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");
    }
}


void DNANMwithRNANMInteraction::check_input_sanity(std::vector<BaseParticle *> &particles){
	DNA2withRNA2Interaction::check_input_sanity(particles);
	//Need to make own function that checks the input sanity
}



void DNANMwithRNANMInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
    int N = particles.size();
    RNANucleotide::set_model(model);

    // checking that the number of nucleotides specified in input file is correct
    if(_nucleotide_types.length() != ndna) {
        throw oxDNAException("Number of nucleotides in the system (%d) doesn't match the length of 'nucleotide_types' (%d)",ndna, _nucleotide_types.length());
    }
    
	if (ndna==0 || ndnas==0) {
        OX_LOG(Logger::LOG_INFO, "No DNA/RNA Particles Specified, Continuing with just Protein Particles");
        for (int i = 0; i < npro; i++) {
            particles[i] = new ACParticle();
        } 
    } else if (npro == 0) {
        OX_LOG(Logger::LOG_INFO, "No Protein Particles Specified, Continuing with just DNA/RNA Particles");
        for (int i = 0; i < ndna; i++) {
            //particles[i] = new DNANucleotide(this->_grooving);
            if(_nucleotide_types[i] == 'D') {
                particles[i] = new DNANucleotide(_grooving_DNA);
            } else {
                particles[i] = new RNANucleotide();
            }   
        }
	} else {
	    if (_firststrand > 0){
            for (int i = 0; i < ndna; i++) {
                //particles[i] = new DNANucleotide(this->_grooving);
                if(_nucleotide_types[i] == 'D') {
                    particles[i] = new DNANucleotide(_grooving_DNA);
                } else {
                    particles[i] = new RNANucleotide();
                }
            }
            for (int i = ndna; i < N; i++) {
                particles[i] = new ACParticle();
            }
	    } else {
            for (int i = 0; i < npro; i++) {
                particles[i] = new ACParticle();
            }
            for (int i = npro; i < N; i++) {
                if(_nucleotide_types[i-npro] == 'D') {
                    particles[i] = new DNANucleotide(_grooving_DNA);    
                } else {
                    particles[i] = new RNANucleotide();
                }
            }
	    }
	}


}


void DNANMwithRNANMInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {

    int N_from_conf = particles.size();
    int my_N, my_N_strands;

    char line[5120];
    std::ifstream topology;
    topology.open(_topology_filename, std::ios::in);

    if (!topology.good())
        throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);

    topology.getline(line, 5120);
    std::stringstream head(line);

    head >> my_N >> my_N_strands >> ndna >> npro >>ndnas;
    if (head.fail()) throw oxDNAException("Problem with header make sure the format is correct for DNANM Interaction");

    if(N_from_conf < 0 || my_N_strands < 0 || my_N_strands > my_N || ndna > my_N || ndna < 0 || npro > my_N || npro < 0 || ndnas < 0 || ndnas > my_N) {
        throw oxDNAException("Problem with header make sure the format is correct for DNANM Interaction");
    }


    int strand, i = 0;
    while (topology.good()) {
        topology.getline(line, 5120);
        if (strlen(line) == 0 || line[0] == '#')
            continue;
        if (i == N_from_conf)
            throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N_from_conf);

        std::stringstream ss(line);
        ss >> strand;
        if (i == 0) {
            _firststrand = strand; //Must be set prior to allocation of particles
            allocate_particles(particles);
            for (int j = 0; j < N_from_conf; j++) {
                particles[j]->index = j;
                particles[j]->type = 0;
                particles[j]->strand_id = 0;
            }
        }

        // Amino Acid
        if (strand < 0) {
            char aminoacid[256];
            int nside, cside;
            ss >> aminoacid >> nside >> cside;

            int x;
            std::set<int> myneighs;
            if (nside >= 0) myneighs.insert(nside);
            if (cside >= 0) myneighs.insert(cside);
            while (ss.good()) {
                ss >> x;
                if (x < 0 || x >= N_from_conf) {
                    throw oxDNAException("Line %d of the topology file has an invalid syntax, neighbor has invalid id",
                                         i + 2);
                }
                myneighs.insert(x);
            }

            ACParticle *p = dynamic_cast< ACParticle * > (particles[i]);

            if (strlen(aminoacid) == 1) {
                p->btype = Utils::decode_aa(aminoacid[0]);
            }

            p->strand_id = abs(strand) + ndnas - 1;
            p->index = i;

            for (std::set<int>::iterator k = myneighs.begin(); k != myneighs.end(); ++k) {
                if (p->index < *k) {
                    p->add_bonded_neighbor(dynamic_cast<ACParticle *> (particles[*k]));
                }
            }
            i++;
        }

        if (strand > 0) {
            char base[256];
            int tmpn3, tmpn5;
            ss >> base >> tmpn3 >> tmpn5;
            
                
            BaseParticle *p = particles[i];   
            /*    
            if(_is_DNA(p)) {
                p = dynamic_cast<DNANucleotide *> (particles[i]);
            } else {
                p = dynamic_cast<RNANucleotide *> (particles[i]);
            }
            */
            
            //DNANucleotide *p = dynamic_cast<DNANucleotide *> (particles[i]);
           


           
            if (tmpn3 < 0) p->n3 = P_VIRTUAL;
            else p->n3 = particles[tmpn3];
            if (tmpn5 < 0) p->n5 = P_VIRTUAL;    
            else p->n5 = particles[tmpn5];

            p->strand_id = strand - 1;
    

            /*
            // allocating particle types for the hybrid interaction
            // -------------------------------------------------------------------------------------------------------------------------
            //making sure that the indexing of 'nucleotide_types' is correct
            if(_nucleotide_types[i-npro] == 'D') {
                p->acid_type = 'D';
            } else {
                p->acid_type = 'R';
            }
            // -------------------------------------------------------------------------------------------------------------------------
            */

            // the base can be either a char or an integer
            if (strlen(base) == 1) {
                p->type = Utils::decode_base(base[0]);
                p->btype = Utils::decode_base(base[0]);

            } else {
                if (atoi(base) > 0) p->type = atoi(base) % 4;
                else p->type = 3 - ((3 - atoi(base)) % 4);
                p->btype = atoi(base);
            }
            


            

            //std::printf("DNA %d %d %d \n", p->index, (p->n3)->index, (p->n5)->index);
            if (p->type == P_INVALID)
                throw oxDNAException("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);

            p->index = i;
            i++;

            // here we fill the affected vector
            if (p->n3 != P_VIRTUAL) p->affected.push_back(ParticlePair(p->n3, p));
            if (p->n5 != P_VIRTUAL) p->affected.push_back(ParticlePair(p, p->n5));
            
            /*
            if(_is_DNA(p)){
                OX_LOG(Logger::LOG_INFO, "D");
            } else {
                OX_LOG(Logger::LOG_INFO, "R");
            }
            */



            //Debug
//            typedef typename std::vector<ParticlePair >::iterator iter;
//            iter it;
//            for (it = p->affected.begin(); it != p->affected.end(); ++it) {
//                printf("Pair %d %d \n", (*it).first->index, (*it).second->index);
//            }

        }
        if (strand == 0) throw oxDNAException("No strand 0 should be present please check topology file");

    }
   
    if (i < N_from_conf)
        throw oxDNAException("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);
    
    topology.close();

    if (my_N != N_from_conf)
        throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");

    *N_strands = my_N_strands; 

   

}



number DNANMwithRNANMInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    if (p->btype >= 0 && q->btype >=0){
        if(p->is_bonded(q)) return this->pair_interaction_bonded(p, q, compute_r, update_forces);
        else return this->pair_interaction_nonbonded(p, q, compute_r, update_forces);
    }
    if ((p->btype >= 0 && q->btype < 0) || (p->btype < 0 && q->btype >= 0)) return this->pair_interaction_nonbonded(p, q, compute_r, update_forces);

    if (p->btype <0 && q->btype <0){
        ACParticle *cp = dynamic_cast< ACParticle * > (p);
        if ((*cp).ACParticle::is_bonded(q)) return this->pair_interaction_bonded(p, q, compute_r, update_forces);
        else return this->pair_interaction_nonbonded(p, q, compute_r, update_forces);
    }
    return 0.f;
}


number DNANMwithRNANMInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(compute_r == true) {
        if (q != P_VIRTUAL && p != P_VIRTUAL) {
            _computed_r = this->_box->min_image(p->pos, q->pos);
        }
    }

    if (p->btype >= 0 && q->btype >=0){
        if(!this->_check_bonded_neighbour(&p, &q, false)) {
            return (number) 0;
        }

        number energy = _na_backbone(p, q, false, update_forces);
        energy += _na_bonded_excluded_volume(p, q, false, update_forces);
        energy += _na_stacking(p, q, false, update_forces);
        return energy;
    }

    if ((p->btype >= 0 && q->btype <0) || (p->btype <0 && q->btype >=0)){
        return 0.f;
    } 

    if ((p->btype <0 && q->btype <0)){
        ACParticle *cp = dynamic_cast< ACParticle * > (p);
        if ((*cp).ACParticle::is_bonded(q)){
            number energy = _protein_spring(p, q, &_computed_r, update_forces);
            energy += _protein_exc_volume(p, q, &_computed_r, update_forces);
            return energy;
        } else{
            return 0.f;
        }
    }

    return 0.f;
}


number DNANMwithRNANMInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if (compute_r == true) {
        _computed_r = this->_box->min_image(p->pos, q->pos);
    }

    if (p->btype >= 0 && q->btype >= 0) { //NA-NA Interaction
        if (_computed_r.norm() >= this->_sqr_rcut) {
            return (number) 0.f;
        } 

        number energy = _na_nonbonded_excluded_volume(p, q, false, update_forces);
        energy += _na_hydrogen_bonding(p, q, false, update_forces);
        energy += _na_cross_stacking(p, q, false, update_forces);
        energy += _na_coaxial_stacking(p, q, false, update_forces);
        energy += _na_debye_huckel(p, q, false, update_forces);
        return energy;
    }

    if ((p->btype >= 0 && q->btype < 0) || (p->btype < 0 && q->btype >= 0)) {
//        if((p->index == 703 && q->index == 482) || (p->index == 482 && q->index == 703)){
//            printf("r.x %.5f r.y %.5f r.z %.5f\n", r->x, r->y, r->z);
//            printf("p.x %.5f p.y %.5f p.z %.5f\n", p->pos.x, p->pos.y, p->pos.z);
//            printf("q.x %.5f q.y %.5f q.z %.5f\n", q->pos.x, q->pos.y, q->pos.z);
//        }
        number energy = _protein_na_exc_volume(p, q, &_computed_r, update_forces);
        return energy;
    }

    if (p->btype < 0 && q->btype < 0) {
        number energy = _protein_exc_volume(p, q, &_computed_r, update_forces);
        return energy;
    }

    return 0.f;
}


number DNANMwithRNANMInteraction::_protein_na_exc_volume(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
  
     BaseParticle *protein;
     BaseParticle *nuc;

     LR_vector force(0, 0, 0);
     LR_vector rcenter = *r;

     if(p->btype >= 0 && q->btype < 0)
     {
         //rcenter = -rcenter;
         protein = q;
         nuc = p;
     }
     else if (p->btype < 0 && q->btype >= 0 )
     {
         rcenter = -rcenter;
         protein = p;
         nuc = q;
     }
     else
     {
         return 0.f;
     }


     LR_vector r_to_back;
     LR_vector r_to_base;
     if(this->_is_DNA(nuc)) {
        r_to_back = rcenter  - nuc->int_centers[DNANucleotide::BACK];
        r_to_base = rcenter  - nuc->int_centers[DNANucleotide::BASE];
     } 
     else {
        r_to_back = rcenter  - nuc->int_centers[RNANucleotide::BACK];
        r_to_base = rcenter  - nuc->int_centers[RNANucleotide::BASE];
     }

     

     LR_vector torquenuc(0,0,0);

     number energy = this->_protein_na_repulsive_lj(r_to_back, force, update_forces, _pro_backbone_sigma, _pro_backbone_b, _pro_backbone_rstar,_pro_backbone_rcut,_pro_backbone_stiffness);
     if (update_forces) {

           

            if(this->_is_DNA(nuc)) {
                torquenuc  -= nuc->int_centers[DNANucleotide::BACK].cross(force);
            } else {
                torquenuc  -= nuc->int_centers[RNANucleotide::BACK].cross(force);
            }
            nuc->force -= force;
            protein->force += force;
     }


     energy += this->_protein_na_repulsive_lj(r_to_base, force, update_forces, _pro_base_sigma, _pro_base_b, _pro_base_rstar, _pro_base_rcut, _pro_base_stiffness);
     if(update_forces) {

         

            if(this->_is_DNA(nuc)) {
                torquenuc  -= nuc->int_centers[DNANucleotide::BASE].cross(force);
            } else {
                torquenuc  -= nuc->int_centers[RNANucleotide::BASE].cross(force);
            }
            nuc->torque += nuc->orientationT * torquenuc;

            nuc->force -= force;
            protein->force += force;
     }


     return energy;
}


number DNANMwithRNANMInteraction::_protein_na_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness) {
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


void DNANMwithRNANMInteraction::init() {
	DNA2withRNA2Interaction::init();
    ndna=0, npro=0, ndnas =0;
	//OLD VERSIONS
    //Backbone-Protein Excluded Volume Parameters
//    _pro_backbone_sigma = 0.4085f;
//    _pro_backbone_rstar= 0.3585f;
//    _pro_backbone_b = 5883.8f;
//    _pro_backbone_rcut = 0.400561f;
//    _pro_backbone_stiffness = 1.0f;
//    //Base-Protein Excluded Volume Parameters
//    _pro_base_sigma = 0.2235f;
//    _pro_base_rstar= 0.1735f;
//    _pro_base_b = 101416.f;
//    _pro_base_rcut = 0.198864f;
//    _pro_base_stiffness = 1.0f;
//    //Protein-Protein Excluded Volume Parameters
//    _pro_sigma = 0.117f;
//    _pro_rstar= 0.087f;
//    _pro_b = 671492.f;
//    _pro_rcut = 0.100161f;
    //NEW_VERSION(Quartic LJ version)
//    _pro_backbone_sigma = 0.68f;
//    _pro_backbone_rstar= 0.679f;
//    _pro_backbone_b = 147802936.f;
//    _pro_backbone_rcut = 0.682945f;
//    _pro_backbone_stiffness = 1.0f;
//    //Base-Protein Excluded Volume Parameters
//    _pro_base_sigma = 0.47f;
//    _pro_base_rstar= 0.45f;
//    _pro_base_b = 157081.f;
//    _pro_base_rcut = 0.506028f;
//    _pro_base_stiffness = 1.0f;
//    //Protein-Protein Excluded Volume Parameters
//    _pro_sigma = 0.55f;
//    _pro_rstar= 0.47f;
//    _pro_b = 80892.1f;
//    _pro_rcut = 0.588787f;
    //that ain't it chief
    //let's try this
    _pro_backbone_sigma = 0.57f;
    _pro_backbone_rstar= 0.569f;
    _pro_backbone_b = 178699253.5f;
    _pro_backbone_rcut = 0.572934f;
    _pro_backbone_stiffness = 1.0f;
    //Base-Protein Excluded Volume Parameters
    _pro_base_sigma = 0.36f;
    _pro_base_rstar= 0.359f;
    _pro_base_b = 296866090.f;
    _pro_base_rcut = 0.362897f;
    _pro_base_stiffness = 1.0f;
    //Protein-Protein Excluded Volume Parameters
    _pro_sigma = 0.35f;
    _pro_rstar = 0.349f;
    _pro_b = 306484596.421f;
    _pro_rcut = 0.352894;

}

//Functions from ACInteraction.h
//Stolen due to inheritance issues



number DNANMwithRNANMInteraction::_protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces) {
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


number DNANMwithRNANMInteraction::_protein_exc_volume(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
    if (p->index != q->index && (p->btype < 0 && q-> btype < 0 )  ){
        LR_vector force(0,0,0);

        number energy =  DNANMwithRNANMInteraction::_protein_repulsive_lj(*r, force, update_forces);

        if(update_forces)
        {
            p->force -= force;
            q->force += force;
        }

        return energy;
    } else {
        return (number) 0.f;
    }
}

number DNANMwithRNANMInteraction::_protein_spring(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
    if(p->btype >= 0 || q->btype >= 0)    //this function is only for proteins
    {
        return 0.f;
    }
    number eqdist;
    char interactiontype;
    std::pair <int, int> keys (std::min(p->index, q->index), std::max(p->index, q->index));
    eqdist = _rknot[keys];
    interactiontype = _potential[keys].first;
    if (eqdist != 0.0) { //only returns number if eqdist is in .par file
        switch (interactiontype) {
            case 's': {
                //Harmonic Spring Potential
                //if ((eqdist < 0.f) || (eqdist > 3.f))  //ensures r0 is nonnegative
//                {
//                    if (keys.first + 1 != keys.second) {
//                        throw oxDNAException("No rknot or invalid rknot value for particle %d and %d rknot was %f",
//                                             q->index, p->index, eqdist);
//                    }
//                }
                number _k = _potential[keys].second; //stiffness of the spring
//                if ((_k == 0) || (_k < 0)) {
//                    throw oxDNAException(
//                            "No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f",
//                            p->index, q->index, _k);
//                }
                number rnorm = r->norm();
                number rinsta = sqrt(rnorm);
                number energy = 0.5 * _k * SQR((rinsta - eqdist));

                if (update_forces) {
                    LR_vector force(*r);
                    force *= (-1.0f * _k) * (rinsta - eqdist) / rinsta;
                    p->force -= force;
                    q->force += force;
//                    printf("p %d q %d d=%f ro=%f df.x=%.8f df.y=%.8f  df.z=%.8f p.x=%.8f p.y=%.8f p.z=%.8f q.x=%.8f "
//                           "q.y=%.8f q.z=%.8f\n", p->index, q->index, rinsta, eqdist, force.x, force.y, force.z, p->force.x,
//                           p->force.y, p->force.z, q->force.x, q->force.y, q->force.z);
                           //"%f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
//                    printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
                }
                return energy;
            }
                break;
            case 'i': {
                //Every Possible Pair of Particles Needs to be Calculated
                number _k = _potential[keys].second; //stiffness of the spring
                if ((_k == 0) || (_k < 0)) {
                    throw oxDNAException(
                            "No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f",
                            p->index, q->index, _k);
                }
                number rnorm = r->norm();
                number rinsta = sqrt(rnorm);
                number energy = 0.5 * _k * SQR(rinsta - eqdist) * (1 / eqdist);

                if (update_forces) {
                    LR_vector force(*r);
                    force *= (-1.0f * _k) * ((rinsta - eqdist) / rinsta) * (1 / eqdist);
                    p->force -= force;
                    q->force += force;
                    //printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
                    //printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
                }
                return energy;
            }
                break;
            case 'e': {
                //Every Possible Pair of Particles Needs to be Calculated
                number _k = _potential[keys].second; //stiffness of the spring
                if ((_k == 0) || (_k < 0)) {
                    throw oxDNAException(
                            "No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f",
                            p->index, q->index, _k);
                }
                number rnorm = r->norm();
                number rinsta = sqrt(rnorm);
                number energy = 0.5 * _k * SQR(rinsta - eqdist) * (1 / (eqdist * eqdist));

                if (update_forces) {
                    LR_vector force(*r);
                    force *= (-1.0f * _k) * ((rinsta - eqdist) / rinsta) * (1 / (eqdist * eqdist));
                    p->force -= force;
                    q->force += force;
                    //printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
                    //printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
                }
                return energy;
            }
                break;
            default: {
                throw oxDNAException(
                        "Interaction type specified in .par file is Invalid, particles %d and %d, switch %c", p->index,
                        q->index, interactiontype);
            }
        }
    } else {
        return (number) 0.f; //returns 0 if no rknot value in parameter value aka they aren't bonded
    }
    //} else {
        //return (number) 0.f; //returns 0 if particle pair consists of particle and itself
    //}

}


number DNANMwithRNANMInteraction::_na_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_backbone(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_bonded_excluded_volume(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_nonbonded_excluded_volume(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_stacking(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_coaxial_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_coaxial_stacking(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_hydrogen_bonding(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_cross_stacking(p, q, false, update_forces);
    } else return 0.f;
}


number DNANMwithRNANMInteraction::_na_debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->btype >= 0 && q->btype >=0){
        return DNA2withRNA2Interaction::_debye_huckel(p, q, false, update_forces);
    } else return 0.f;
}











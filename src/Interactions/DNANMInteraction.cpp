/*
 * DNANMInteraction.cpp
 *
 *  Created on: Feb 2026
 *  Author: jonah,subhajit
 *  Inherits from the DNA2Interaction Class
 *  Uses DNA2 Model and ANM Protein Model
 */


#include "DNANMInteraction.h"

#include "../Particles/DNANucleotide.h"
#include "../Particles/ANMParticle.h"
#include "../Utilities/Utils.h"

#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <map>
#include <set>
#include <utility>
#include <cstring>


namespace {

std::string trim_copy(const std::string &value) {
    std::size_t first = 0;
    while(first < value.size() &&
          std::isspace(static_cast<unsigned char>(value[first]))) {
        first++;
    }

    std::size_t last = value.size();
    while(last > first &&
          std::isspace(static_cast<unsigned char>(value[last - 1]))) {
        last--;
    }

    return value.substr(first, last - first);
}

bool starts_with(const std::string &value, const std::string &prefix) {
    return value.size() >= prefix.size() &&
           value.compare(0, prefix.size(), prefix) == 0;
}

double parse_finite_double_text(
    const std::string &value,
    const std::string &what
) {
    if(value.empty()) {
        throw oxDNAException(
            "%s is empty",
            what.c_str()
        );
    }

    // oxDNA Release builds use -ffast-math, which includes
    // -ffinite-math-only. Reject NaN/Inf spellings lexically before
    // conversion rather than relying on floating-point classification.
    std::string lower_value(value);

    std::transform(
        lower_value.begin(),
        lower_value.end(),
        lower_value.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        }
    );

    std::size_t token_start = 0;

    if(!lower_value.empty() &&
       (lower_value[0] == '+' || lower_value[0] == '-')) {
        token_start = 1;
    }

    const std::string unsigned_value =
        lower_value.substr(token_start);

    if(unsigned_value == "inf" ||
       unsigned_value == "infinity" ||
       starts_with(unsigned_value, "nan")) {
        throw oxDNAException(
            "%s is out of range or non-finite (`%s')",
            what.c_str(),
            value.c_str()
        );
    }

    errno = 0;
    char *end = nullptr;
    const char *begin = value.c_str();
    const double parsed = std::strtod(begin, &end);

    if(end == begin) {
        throw oxDNAException(
            "%s has invalid value `%s'",
            what.c_str(),
            value.c_str()
        );
    }

    while(end != nullptr &&
          *end != '\0' &&
          std::isspace(static_cast<unsigned char>(*end))) {
        end++;
    }

    if(end == nullptr || *end != '\0') {
        throw oxDNAException(
            "%s has invalid trailing characters in value `%s'",
            what.c_str(),
            value.c_str()
        );
    }

    if(errno == ERANGE) {
        throw oxDNAException(
            "%s is out of range or non-finite (`%s')",
            what.c_str(),
            value.c_str()
        );
    }

    return parsed;
}

number checked_number_product(
    number lhs,
    number rhs,
    const char *what
) {
    const number lhs_abs = std::fabs(lhs);
    const number rhs_abs = std::fabs(rhs);
    const number limit = std::numeric_limits<number>::max();

    if(lhs_abs != (number) 0.0 &&
       rhs_abs > limit / lhs_abs) {
        throw oxDNAException(
            "%s overflows the compiled CPU numerical type",
            what
        );
    }

    return lhs * rhs;
}

number checked_reduced_protein_charge(
    number protein_charge_scale,
    double charge_e,
    number dna_effective_charge_e,
    int local_index
) {
    const double scale = static_cast<double>(protein_charge_scale);
    const double dna_charge = static_cast<double>(dna_effective_charge_e);
    const double product_limit = std::numeric_limits<double>::max();

    const double scale_abs = std::fabs(scale);
    const double charge_abs = std::fabs(charge_e);

    if(scale_abs != 0.0 &&
       charge_abs > product_limit / scale_abs) {
        throw oxDNAException(
            "Reduced protein charge for local index %d overflows while multiplying protein_charge_scale by charge_e",
            local_index
        );
    }

    const double numerator = scale * charge_e;
    const double numerator_abs = std::fabs(numerator);

    // dna_effective_charge_e is required to be positive. Division can
    // overflow only when the divisor is below one; the product on the
    // right remains <= DBL_MAX and therefore cannot overflow here.
    if(dna_charge < 1.0 &&
       numerator_abs > product_limit * dna_charge) {
        throw oxDNAException(
            "Reduced protein charge for local index %d overflows while dividing by dna_effective_charge_e",
            local_index
        );
    }

    const double reduced = numerator / dna_charge;
    const double number_limit =
        static_cast<double>(std::numeric_limits<number>::max());

    if(reduced < -number_limit || reduced > number_limit) {
        throw oxDNAException(
            "Reduced protein charge for local index %d is outside the compiled CPU numerical range",
            local_index
        );
    }

    return static_cast<number>(reduced);
}

bool parse_optional_number_strict(
    input_file &inp,
    const char *key,
    number &dest
) {
    std::string raw;

    if(getInputString(&inp, key, raw, 0) != KEY_FOUND) {
        return false;
    }

    const std::string value = trim_copy(raw);
    const std::string what =
        std::string("Numeric input key `") + key + "'";

    const double parsed =
        parse_finite_double_text(value, what);

    const double limit =
        static_cast<double>(std::numeric_limits<number>::max());

    if(parsed < -limit || parsed > limit) {
        throw oxDNAException(
            "Numeric input key `%s' is outside the compiled CPU numerical range (`%s')",
            key,
            value.c_str()
        );
    }

    dest = static_cast<number>(parsed);
    return true;
}

} // namespace

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

    _protein_dna_electrostatics = false;
    _protein_protein_electrostatics = false;
    _exclude_intragroup_electrostatics = true;

    _protein_charge_scale = (number) 1.0;
    _protein_dna_dh_scale = (number) 1.0;
    _protein_protein_dh_scale = (number) 1.0;

    // Standard oxDNA2 effective backbone charge reference.
    _dna_effective_charge_e = (number) 0.815;

    _protein_charge_format_version = 0;


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

    getInputBool(
        &inp,
        "protein_dna_electrostatics",
        &_protein_dna_electrostatics,
        0
    );



    getInputBool(
        &inp,
        "exclude_intragroup_electrostatics",
        &_exclude_intragroup_electrostatics,
        0
    );

    parse_optional_number_strict(
        inp,
        "protein_charge_scale",
        _protein_charge_scale
    );

    parse_optional_number_strict(
        inp,
        "protein_dna_dh_scale",
        _protein_dna_dh_scale
    );

    parse_optional_number_strict(
        inp,
        "protein_protein_dh_scale",
        _protein_protein_dh_scale
    );

    parse_optional_number_strict(
        inp,
        "dna_effective_charge_e",
        _dna_effective_charge_e
    );

    if(_dna_effective_charge_e <= (number) 0.0) {
        throw oxDNAException(
            "dna_effective_charge_e must be positive"
        );
    }

    if(_protein_charge_scale < (number) 0.0) {
        throw oxDNAException(
            "protein_charge_scale must be non-negative"
        );
    }

    if(_protein_dna_dh_scale < (number) 0.0) {
        throw oxDNAException(
            "protein_dna_dh_scale must be non-negative"
        );
    }

    if(_protein_protein_dh_scale < (number) 0.0) {
        throw oxDNAException(
            "protein_protein_dh_scale must be non-negative"
        );
    }

    if(_protein_dna_electrostatics) {
        ADD_INTERACTION_TO_MAP(
            PRO_DNA_DH,
            _protein_dna_debye_huckel
        );
    }



    getInputBool(&inp, "protein_protein_electrostatics", &_protein_protein_electrostatics, 0);
    if(_protein_protein_electrostatics) {
        ADD_INTERACTION_TO_MAP(PRO_PRO_DH, _protein_protein_debye_huckel);
    }
    if(_protein_dna_electrostatics ||
       _protein_protein_electrostatics) {

        if(getInputString(
               &inp,
               "protein_charge_file",
               _protein_charge_file,
               0
           ) != KEY_FOUND) {

            throw oxDNAException(
                "protein_charge_file is required when protein electrostatics are enabled"
            );
        }

        OX_LOG(
            Logger::LOG_INFO,
            "Protein electrostatics enabled: PD=%d PP=%d, "
            "charge scale=%g, PD scale=%g, PP scale=%g, "
            "oxDNA2 charge-normalisation reference=%g e",
            (int) _protein_dna_electrostatics,
            (int) _protein_protein_electrostatics,
            _protein_charge_scale,
            _protein_dna_dh_scale,
            _protein_protein_dh_scale,
            _dna_effective_charge_e
        );
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
                auto *q = dynamic_cast<ANMParticle *>(particles[key1]);
                if (key1 < key2) {
                    q->add_bonded_neighbor(particles[key2]);
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
        for (int i = 0; i < npro; i++) particles[i] = new ANMParticle(_angular);

    } else if (npro == 0) {
        OX_LOG(Logger::LOG_INFO, "No Protein Particles Specified, Continuing with just DNA Particles");
        for (int i = 0; i < ndna; i++) particles[i] = new DNANucleotide(this->_grooving);

    } else {
        if (_firststrand > 0){
            for (int i = 0; i < ndna; i++) particles[i] = new DNANucleotide(this->_grooving);
            // Protein
            for (uint i = ndna; i < particles.size(); i++) particles[i] = new ANMParticle(_angular);
        } else {
            // Protein
            for (int i = 0; i < npro; i++) particles[i] = new ANMParticle(_angular);

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
    if (head.fail()) throw oxDNAException("Invalid topology header for DNANM. This interaction currently supports the mixed legacy format: '<N> <N_strands> <ndna> <npro> <ndnas>'.");

    if(my_N_strands < 0 || my_N_strands > my_N || ndna > my_N || ndna < 0 || npro > my_N || npro < 0 || ndnas < 0 || ndnas > my_N) {
        throw oxDNAException("Invalid topology header for DNANM. This interaction currently supports the mixed legacy format: '<N> <N_strands> <ndna> <npro> <ndnas>'.");
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
            auto *p = static_cast<ANMParticle*>(particles[i]); //Subho
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
            p->invmass = 1.f/p->mass;

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

            // p->mass = 1.f; //subho not needed as DNA is 1 by default
            // p->invmass = 1.f;

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
    if(_protein_dna_electrostatics || _protein_protein_electrostatics) {
        _prepare_protein_charge_topology_metadata(my_N);
        _load_protein_charge_file(my_N);
    }
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
        if (rnorm >= (_protein_protein_electrostatics ? _pro_pro_sqr_rcut : _pro_sqr_rcut)) return (number) 0.f;
        number energy = _protein_exc_volume(p, q, compute_r, update_forces);
        if(_protein_protein_electrostatics)
            energy += _protein_protein_debye_huckel(p, q, false, update_forces);
        return energy;
//        return 0.f;
    }else if(interaction_type == 1) { //protein-dna
        if (rnorm >= _pro_dna_sqr_rcut) return (number) 0.f;
        number energy = _protein_dna_exc_volume(p, q, compute_r, update_forces);
        if(_protein_dna_electrostatics)
            energy += _protein_dna_debye_huckel(p, q, false, update_forces);
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
    _pro_pro_sqr_rcut = _pro_sqr_rcut;
    _update_electrostatic_cutoffs();
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







// Optional, strictly opt-in protein screened electrostatics.

void DNANMInteraction::_prepare_protein_charge_topology_metadata(int total_particles) {
    _protein_local_to_global.clear();
    _protein_topology_strand.clear();
    _protein_residue_code.clear();
    std::ifstream in(_topology_filename);
    std::string row;
    if(!in || !std::getline(in, row))
        throw oxDNAException("Cannot read protein topology for electrostatics");
    int n = 0, nstrands = 0, dna = 0, pro = 0, dnas = 0;
    std::istringstream header(row);
    std::string extra;
    if(!(header >> n >> nstrands >> dna >> pro >> dnas) ||
       (header >> extra) || n != total_particles || pro != npro || dna != ndna ||
       n != dna + pro || pro <= 0)
        throw oxDNAException("Invalid protein electrostatics topology header");
    const bool protein_first = (_firststrand < 0);
    int index = 0;
    while(std::getline(in, row)) {
        const std::string line=trim_copy(row);
        if(line.empty() || line[0]=='#') continue;
        if(index >= n) throw oxDNAException("Extra topology row in protein charge mapping");
        std::istringstream fields(line);
        int strand=0;
        if(!(fields >> strand) || strand==0)
            throw oxDNAException("Invalid strand in protein charge mapping at index %d",index);
        const bool expected_protein=protein_first ? index<pro : index>=dna;
        if((strand<0)!=expected_protein)
            throw oxDNAException("Interleaved DNA/protein topology unsupported in protein electrostatics");
        if(expected_protein) {
            std::string residue;
            int n3, n5;
            if(!(fields >> residue >> n3 >> n5) || residue.size()!=1)
                throw oxDNAException("Invalid protein topology row at index %d",index);
            _protein_local_to_global.push_back(index);
            _protein_topology_strand.push_back(strand);
            _protein_residue_code.push_back(static_cast<char>(
                std::toupper(static_cast<unsigned char>(residue[0]))));
        }
        ++index;
    }
    if(index!=n || static_cast<int>(_protein_local_to_global.size())!=pro)
        throw oxDNAException("Incorrect number of protein topology rows");
}


void DNANMInteraction::_load_protein_charge_file(int total_particles) {
    if(!_protein_dna_electrostatics && !_protein_protein_electrostatics) return;
    if(npro <= 0 || static_cast<int>(_protein_local_to_global.size()) != npro)
        throw oxDNAException("Invalid protein electrostatic topology mapping");
    std::ifstream in(_protein_charge_file);
    if(!in) throw oxDNAException("Cannot open protein charge file '%s'", _protein_charge_file.c_str());
    struct Entry { int local, strand, group; char residue; double charge; };
    std::vector<Entry> entries;
    std::map<std::string, std::string> meta;
    std::string line;
    bool header_seen=false, saw_data=false;
    int version=0, line_no=0;
    while(std::getline(in,line)) {
        ++line_no;
        line=trim_copy(line);
        if(line.empty()) continue;
        if(line[0]=='#') {
            const std::string comment=trim_copy(line.substr(1));
            if(starts_with(comment,"DNANM_ES_CHARGES")) {
                if(header_seen || saw_data)
                    throw oxDNAException("Duplicate or misplaced charge-format header on line %d",line_no);
                std::istringstream header(comment);
                std::string magic, trailing;
                if(!(header>>magic>>version) || magic!="DNANM_ES_CHARGES" ||
                   (header>>trailing) || version!=1)
                    throw oxDNAException("Invalid versioned charge header on line %d",line_no);
                header_seen=true;
            } else {
                const auto eq=comment.find('=');
                if(eq!=std::string::npos) {
                    const std::string key=trim_copy(comment.substr(0,eq));
                    const std::string value=trim_copy(comment.substr(eq+1));
                    if(!meta.emplace(key,value).second)
                        throw oxDNAException("Duplicate charge metadata '%s'",key.c_str());
                }
            }
            continue;
        }
        saw_data=true;
        const auto hash=line.find('#');
        if(hash!=std::string::npos) line=trim_copy(line.substr(0,hash));
        if(line.empty()) continue;
        Entry r={-1,0,-1,'\0',0.0};
        std::string residue, charge_text, extra;
        std::istringstream record(line);
        if(version==1) {
            if(!(record>>r.local>>r.strand>>residue>>charge_text>>r.group) ||
               residue.size()!=1 || (record>>extra))
                throw oxDNAException("Malformed version-1 charge record on line %d",line_no);
            r.residue=static_cast<char>(std::toupper(static_cast<unsigned char>(residue[0])));
        } else {
            if(!(record>>r.local>>charge_text>>r.group) || (record>>extra))
                throw oxDNAException("Malformed legacy charge record on line %d",line_no);
        }
        r.charge=parse_finite_double_text(charge_text,
                  "protein charge line "+std::to_string(line_no));
        entries.push_back(r);
    }
    if(in.bad()) throw oxDNAException("I/O error in protein charge file");
    auto require_meta=[&](const char *key)->std::string {
        auto it=meta.find(key);
        if(it==meta.end()) throw oxDNAException("Charge file missing metadata '%s'",key);
        return it->second;
    };
    auto parse_meta_int=[&](const std::string &v, const char *key)->int {
        std::istringstream ss(v);
        int n; std::string trailing;
        if(!(ss>>n) || (ss>>trailing))
            throw oxDNAException("Invalid integer metadata '%s'",key);
        return n;
    };
    if(version==1) {
        if(require_meta("index_space")!="protein_local" || require_meta("index_base")!="0")
            throw oxDNAException("Charge indices must be zero-based protein_local");
        if(parse_meta_int(require_meta("n_protein_beads"),"n_protein_beads")!=npro)
            throw oxDNAException("Charge file topology protein count mismatch");
        _protein_charge_model=meta.count("charge_model") ? meta.at("charge_model") : "unspecified";
        _protein_charge_source=meta.count("charge_source") ? meta.at("charge_source") : "unspecified";
    } else {
        _protein_charge_model="legacy-unspecified";
        _protein_charge_source="legacy-unspecified";
        OX_LOG(Logger::LOG_WARNING,"Legacy three-column protein charge format used");
    }
    if(static_cast<int>(entries.size())!=npro)
        throw oxDNAException("Charge file has %d entries, expected %d",(int)entries.size(),npro);
    std::vector<char> seen(npro,0);
    _protein_charge_reduced.assign(total_particles,(number)0);
    _protein_group.assign(total_particles,-1);
    std::set<int> groups;
    for(const auto &e:entries) {
        if(e.local<0 || e.local>=npro || seen[e.local])
            throw oxDNAException("Duplicate or invalid protein charge local index %d",e.local);
        if(e.group<0) throw oxDNAException("Negative charge group for protein %d",e.local);
        if(version==1 && (e.strand!=_protein_topology_strand[e.local] ||
                         e.residue!=_protein_residue_code[e.local]))
            throw oxDNAException("Protein charge/topology strand or residue mismatch at %d",e.local);
        const int global=_protein_local_to_global[e.local];
        if(global<0 || global>=total_particles)
            throw oxDNAException("Invalid global charge mapping at %d",e.local);
        _protein_charge_reduced[global]=checked_reduced_protein_charge(
            _protein_charge_scale,e.charge,_dna_effective_charge_e,e.local);
        _protein_group[global]=e.group;
        seen[e.local]=1;
        groups.insert(e.group);
    }
    if(version==1 && meta.count("n_groups") &&
       parse_meta_int(meta.at("n_groups"),"n_groups")!=static_cast<int>(groups.size()))
        throw oxDNAException("Charge metadata n_groups does not match file");
    _protein_charge_format_version=version;
    OX_LOG(Logger::LOG_INFO,"Loaded %d protein charges from %s (format=%d)",
           npro,_protein_charge_file.c_str(),version);
}

void DNANMInteraction::_update_electrostatic_cutoffs() {

    _pro_pro_sqr_rcut =
        _pro_sqr_rcut;

    number required_sqr_cut =
        (number) 0.0;

    if(_protein_protein_electrostatics &&
       _protein_protein_dh_scale >
           (number) 0.0) {

        _pro_pro_sqr_rcut =
            std::max(
                _pro_pro_sqr_rcut,
                SQR(_debye_huckel_RC)
            );

        required_sqr_cut =
            std::max(
                required_sqr_cut,
                _pro_pro_sqr_rcut
            );
    }

    if(_protein_dna_electrostatics &&
       _protein_dna_dh_scale >
           (number) 0.0) {

        const number backbone_offset =
            (_grooving)
            ? sqrt(
                SQR(POS_MM_BACK1)
                +
                SQR(POS_MM_BACK2)
              )
            : std::fabs(
                (number) POS_BACK
              );

        const number pd_com_cut =
            backbone_offset
            +
            _debye_huckel_RC;

        _pro_dna_sqr_rcut =
            std::max(
                _pro_dna_sqr_rcut,
                SQR(pd_com_cut)
            );

        required_sqr_cut =
            std::max(
                required_sqr_cut,
                _pro_dna_sqr_rcut
            );
    }

    if(required_sqr_cut >
       (number) 0.0) {

        const number required_cut =
            sqrt(required_sqr_cut);

        if(required_cut > _rcut) {

            _rcut =
                required_cut;

            _sqr_rcut =
                SQR(_rcut);
        }
    }
}

number DNANMInteraction::_scaled_debye_huckel(
    const LR_vector &r,
    number charge_product,
    LR_vector &force,
    bool update_forces
) {

    force = LR_vector(0, 0, 0);

    if(charge_product ==
       (number) 0.0) {
        return (number) 0.0;
    }

    const number rmod =
        r.module();

    if(rmod <=
           std::numeric_limits<number>::epsilon()
       ||
       rmod >=
           _debye_huckel_RC) {

        return (number) 0.0;
    }

    number energy;

    if(rmod <
       _debye_huckel_RHIGH) {

        const number exponential =
            exp(
                _minus_kappa *
                rmod
            );

        energy =
            charge_product
            *
            exponential
            *
            (_debye_huckel_prefactor /
             rmod);

        if(update_forces) {

            force =
                (r / rmod)
                *
                (
                    charge_product
                    *
                    (
                        -_debye_huckel_prefactor
                        *
                        exponential
                    )
                    *
                    (
                        _minus_kappa /
                        rmod
                        -
                        (number) 1.0 /
                        SQR(rmod)
                    )
                );
        }
    }
    else {

        energy =
            charge_product
            *
            _debye_huckel_B
            *
            SQR(
                rmod -
                _debye_huckel_RC
            );

        if(update_forces) {

            force =
                -(r / rmod)
                *
                (
                    charge_product
                    *
                    (number) 2.0
                    *
                    _debye_huckel_B
                    *
                    (
                        rmod -
                        _debye_huckel_RC
                    )
                );
        }
    }

    return energy;
}

number DNANMInteraction::_protein_dna_debye_huckel(
    BaseParticle *p,
    BaseParticle *q,
    bool compute_r,
    bool update_forces
) {

    if(!_protein_dna_electrostatics) {
        return (number) 0.0;
    }

    BaseParticle *protein = nullptr;
    BaseParticle *nuc = nullptr;

    LR_vector rcenter =
        _computed_r;

    const int pid =
        get_id(p->btype);

    const int qid =
        get_id(q->btype);

    if(pid == 0 && qid == 1) {

        nuc = p;
        protein = q;
    }
    else if(pid == 1 &&
            qid == 0) {

        protein = p;
        nuc = q;

        rcenter =
            -rcenter;
    }
    else {
        return (number) 0.0;
    }

    if(protein->index < 0 ||
       protein->index >=
           (int) _protein_charge_reduced.size()) {

        throw oxDNAException(
            "DNANM protein-DNA electrostatics requested for protein particle %d outside charge array of size %d",
            protein->index,
            (int) _protein_charge_reduced.size()
        );
    }

    // DNA backbone carries negative effective charge.
    number charge_product =
        checked_number_product(
            -_protein_dna_dh_scale,
            _protein_charge_reduced[
                protein->index
            ],
            "Protein-DNA electrostatic charge factor"
        );

    if(_debye_huckel_half_charged_ends &&
       (
           nuc->n3 == P_VIRTUAL ||
           nuc->n5 == P_VIRTUAL
       )) {

        charge_product *=
            (number) 0.5;
    }

    const LR_vector rback =
        rcenter
        -
        nuc->int_centers[
            DNANucleotide::BACK
        ];

    LR_vector force_on_protein(
        0,
        0,
        0
    );

    const number energy =
        _scaled_debye_huckel(
            rback,
            charge_product,
            force_on_protein,
            update_forces
        );

    if(update_forces &&
       energy != (number) 0.0) {

        const LR_vector force_on_nuc =
            -force_on_protein;

        nuc->force +=
            force_on_nuc;

        protein->force +=
            force_on_protein;

        const LR_vector torque_nuc =
            nuc->int_centers[
                DNANucleotide::BACK
            ].cross(
                force_on_nuc
            );

        nuc->torque +=
            nuc->orientationT
            *
            torque_nuc;

        // Current upstream particle-aware stress API.
        // rcenter = protein COM - nucleotide COM.
        // force_on_protein = force acting on second particle.
        _update_stress_tensor(
            nuc,
            protein,
            rcenter,
            force_on_protein
        );
    }

    return energy;
}

number DNANMInteraction::_protein_protein_debye_huckel(
    BaseParticle *p,
    BaseParticle *q,
    bool compute_r,
    bool update_forces
) {

    if(!_protein_protein_electrostatics ||
       get_id(p->btype) != 1 ||
       get_id(q->btype) != 1) {

        return (number) 0.0;
    }

    // Direct term API and nonbonded lists must agree: bonded ANM springs
    // contribute through bonded mechanics, never screened PP DH.
    if(p->is_bonded(q)) return (number) 0.0;
    if(p->index < 0 ||
       q->index < 0 ||
       p->index >=
           (int) _protein_charge_reduced.size() ||
       q->index >=
           (int) _protein_charge_reduced.size()) {

        throw oxDNAException(
            "DNANM protein-protein electrostatics requested for particles %d/%d outside charge array of size %d",
            p->index,
            q->index,
            (int) _protein_charge_reduced.size()
        );
    }

    if(_exclude_intragroup_electrostatics &&
       _protein_group[p->index] ==
       _protein_group[q->index]) {

        return (number) 0.0;
    }

    const number scaled_p_charge =
        checked_number_product(
            _protein_protein_dh_scale,
            _protein_charge_reduced[p->index],
            "Protein-protein electrostatic scale/charge product"
        );

    const number charge_product =
        checked_number_product(
            scaled_p_charge,
            _protein_charge_reduced[q->index],
            "Protein-protein electrostatic charge product"
        );

    LR_vector force_on_q(
        0,
        0,
        0
    );

    const number energy =
        _scaled_debye_huckel(
            _computed_r,
            charge_product,
            force_on_q,
            update_forces
        );

    if(update_forces &&
       energy != (number) 0.0) {

        p->force -=
            force_on_q;

        q->force +=
            force_on_q;

        // Current upstream particle-aware stress API.
        _update_stress_tensor(
            p,
            q,
            _computed_r,
            force_on_q
        );
    }

    return energy;
}

number DNANMInteraction::pair_interaction_term(
    int name,
    BaseParticle *p,
    BaseParticle *q,
    bool compute_r,
    bool update_forces
) {

    const int interaction_type =
        get_id(p->btype)
        +
        get_id(q->btype);

    const bool bonded =
        p->is_bonded(q);

    BaseParticle *left = p;
    BaseParticle *right = q;

    // A supplied cache is oriented from the caller's p to q. Canonical PP
    // endpoints require the opposite vector when that order is reversed.
    // Restore the caller's cache on every return and on exceptions so repeated
    // term queries cannot alternately flip or silently invalidate its geometry.
    const bool reverse_cached_r = interaction_type == 2 &&
        left->index > right->index && !compute_r;
    struct CachedDisplacementGuard {
        LR_vector &target;
        LR_vector saved;
        bool restore;
        ~CachedDisplacementGuard() {
            if(restore) target = saved;
        }
    } cached_displacement_guard{
        _computed_r,
        reverse_cached_r ? _computed_r : LR_vector(0, 0, 0),
        reverse_cached_r
    };

    if(interaction_type == 2 &&
       left->index > right->index) {

        std::swap(left, right);
        if(reverse_cached_r) _computed_r = -_computed_r;
    }

    if(compute_r) {
        _computed_r =
            this->_box->min_image(
                left->pos,
                right->pos
            );
    }

    if(interaction_type == 0) {

        const bool bonded_term =
            name == BACKBONE
            ||
            name == BONDED_EXCLUDED_VOLUME
            ||
            name == STACKING;

        const bool nonbonded_term =
            name == NONBONDED_EXCLUDED_VOLUME
            ||
            name == HYDROGEN_BONDING
            ||
            name == CROSS_STACKING
            ||
            name == COAXIAL_STACKING
            ||
            name == DEBYE_HUCKEL;

        if(
            (bonded && bonded_term)
            ||
            (!bonded && nonbonded_term)
        ) {

            return BaseInteraction::
                pair_interaction_term(
                    name,
                    p,
                    q,
                    false,
                    update_forces
                );
        }

        return (number) 0.0;
    }

    if(interaction_type == 1) {

        if(!bonded &&
           (
               name == PRO_DNA_EXC_VOL
               ||
               (
                   name == PRO_DNA_DH
                   &&
                   _protein_dna_electrostatics
               )
           )) {

            return BaseInteraction::
                pair_interaction_term(
                    name,
                    p,
                    q,
                    false,
                    update_forces
                );
        }

        return (number) 0.0;
    }

    if(interaction_type == 2) {

        if(name == SPRING) {
            return bonded
                ? _protein_spring(
                    left,
                    right,
                    false,
                    update_forces
                  )
                : (number) 0.0;
        }

        if(name == PRO_EXC_VOL) {
            return _protein_exc_volume(
                left,
                right,
                false,
                update_forces
            );
        }

        if(name == PRO_PRO_DH) {

            if(bonded ||
               !_protein_protein_electrostatics) {

                return (number) 0.0;
            }

            return _protein_protein_debye_huckel(
                left,
                right,
                false,
                update_forces
            );
        }

        if(name == PRO_ANG_POT) {

            if(!_angular ||
               !bonded ||
               right->index -
               left->index != 1) {

                return (number) 0.0;
            }

            return _protein_ang_pot(
                left,
                right,
                false,
                update_forces
            );
        }
    }

    return (number) 0.0;
}

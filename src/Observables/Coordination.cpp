/*
 * Coordination.cpp
 *
 * Created on: 10/17/2025
 *     Author: Lorenzo
*/

#include "Coordination.h"

Coordination::Coordination() {

}

Coordination::~Coordination() {

}

void Coordination::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

    getInputString(&my_inp, "op_file", _op_file, 1);

    if(getInputNumber(&my_inp, "d0", &_d0, 0) == KEY_NOT_FOUND) {
        _d0 = 1.2;
    }
	if(getInputNumber(&my_inp, "r0", &_r0, 0) == KEY_NOT_FOUND) {
        _r0 = 1.2;
    }
    if(getInputInt(&my_inp, "n", &_n, 0) == KEY_NOT_FOUND) {
        _n = 6;
    }
}

void Coordination::init() {
	BaseObservable::init();

    if(_n % 2 != 0) {
        throw oxDNAException("Coordination: exponent n must be an even integer");
    }

    OrderParameters op;
    op.init_from_file(_op_file.c_str(), CONFIG_INFO->particles(), CONFIG_INFO->N());
    if(op.get_distance_parameters_count() > 0) {
        throw oxDNAException("Coordination: distance-based order parameters are not supported");
    }

    if(op.get_hb_parameters_count() == 0) {
        throw oxDNAException("Coordination: no hydrogen-bond-based order parameters found in the specified op_file");
    }

    if(op.get_hb_parameters_count() > 1) {
        throw oxDNAException("Coordination: only one hydrogen-bond-based order parameter is supported");
    }

    auto all_index_pairs = op.get_hb_particle_list();
    std::vector<int> p_indices;
    p_indices.reserve(all_index_pairs.size() * 2);
    for(auto &pair : all_index_pairs) {
        p_indices.push_back(pair.first);
        p_indices.push_back(pair.second);

        BaseParticle *p1 = CONFIG_INFO->particles()[pair.first];
        BaseParticle *p2 = CONFIG_INFO->particles()[pair.second];

        _all_pairs.push_back(std::make_pair(p1, p2));
    }

    // sorting is not necessary, but nice if we want to print all indices for debugging purposes
    // and also convenient to check for duplicates
    std::sort(p_indices.begin(), p_indices.end());

    // here we look for duplicates
    auto it = std::adjacent_find(p_indices.begin(), p_indices.end());
    if(it != p_indices.end()) {
        throw oxDNAException("LTCoordination: duplicate particle index found in the op_file: %d", *it);
    }
}

std::string Coordination::get_output_string(llint curr_step) {
    number coordination = 0.0;
    for(auto &pair : _all_pairs) {
        number r = std::sqrt(CONFIG_INFO->box->sqr_min_image_distance(pair.first->pos, pair.second->pos));
        coordination += 1.0 / (1.0 + std::pow((r - _d0) / _r0, _n));
    }

    return Utils::sformat("%lf", coordination);
}
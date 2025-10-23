/*
 * LTCoordination.cpp
 *
 * Created on: 10/16/2025
 *     Author: Lorenzo
*/

#include "LTCoordination.h"

#include "meta_utils.h"
#include "../../Utilities/OrderParameters.h"
#include "../../Utilities/Utils.h"

LTCoordination::LTCoordination() {

}

LTCoordination::~LTCoordination() {

}

std::tuple<std::vector<int>, std::string> LTCoordination::init(input_file &inp) {
    BaseForce::init(inp);

    std::string op_file;
    getInputString(&inp, "op_file", op_file, 1);

    OrderParameters op;
    op.init_from_file(op_file.c_str(), CONFIG_INFO->particles(), CONFIG_INFO->N());

    if(op.get_distance_parameters_count() > 0) {
        throw oxDNAException("LTCoordination: distance-based order parameters are not supported");
    }

    if(op.get_hb_parameters_count() == 0) {
        throw oxDNAException("LTCoordination: no hydrogen-bond-based order parameters found in the specified op_file");
    }

    if(op.get_hb_parameters_count() > 1) {
        throw oxDNAException("LTCoordination: only one hydrogen-bond-based order parameter is supported");
    }

    auto all_index_pairs = op.get_hb_particle_list();
    std::vector<int> p_indices;
    p_indices.reserve(all_index_pairs.size() * 2);
    for(auto &pair : all_index_pairs) {
        p_indices.push_back(pair.first);
        p_indices.push_back(pair.second);

        BaseParticle *p1 = CONFIG_INFO->particles()[pair.first];
        BaseParticle *p2 = CONFIG_INFO->particles()[pair.second];
        particle_to_pair_index[p1] = all_pairs.size();
        particle_to_pair_index[p2] = all_pairs.size();

        all_pairs.push_back(std::make_pair(p1, p2));
    }

    // sorting is not necessary, but nice if we want to print all indices for debugging purposes
    // and also convenient to check for duplicates
    std::sort(p_indices.begin(), p_indices.end());

    // here we look for duplicates
    auto it = std::adjacent_find(p_indices.begin(), p_indices.end());
    if(it != p_indices.end()) {
        throw oxDNAException("LTCoordination: duplicate particle index found in the op_file: %d", *it);
    }

    getInputNumber(&inp, "d0", &d0, 0);
	getInputNumber(&inp, "r0", &r0, 0);
    getInputInt(&inp, "n", &n, 0);

    if(n % 2 != 0) {
        throw oxDNAException("LTCoordination: exponent n must be an even integer");
    }

    // build the grid
    if(getInputNumber(&inp, "coord_min", &coord_min, 0) == KEY_NOT_FOUND) {
        coord_min = 0.0;
    }
	if(getInputNumber(&inp, "coord_max", &coord_max, 0) == KEY_NOT_FOUND) {
        coord_max = all_index_pairs.size() * 1.01; // maximum coordination number possible + small buffer
    }

	getInputInt(&inp, "N_grid", &N_grid, 1);
    d_coord = (coord_max - coord_min) / (N_grid - 1.0);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	potential_grid = meta::split_to_numbers(potential_string, ",");

    if(potential_grid.size() != (size_t) N_grid) {
        throw oxDNAException("LTCoordination: potential_grid size (%u) != N_grid (%d)", potential_grid.size(), N_grid);
    }

    return std::make_tuple(p_indices, "LTCoordination force");
}

LR_vector LTCoordination::value(llint step, LR_vector &pos) {
    auto pair = all_pairs[particle_to_pair_index[_current_particle]];
    int sign = (pair.first == _current_particle) ? -1 : 1;

    LR_vector r_vec = CONFIG_INFO->box->min_image(pair.first->pos, pair.second->pos);
    number r = r_vec.norm();
    double x = (r - d0) / r0;
    double xn = std::pow(x, n);
    double dcoord_dr = -(n / r0) * std::pow(x, n - 1) / (SQR(1.0 + xn));

    number coord = Utils::clamp(_coordination(), coord_min, coord_max); // ensure we are within the grid limits
    int ic_left = std::floor((coord - coord_min) / d_coord);
	int ic_right = ic_left + 1;

    number df_dcoord = 0;
	if((ic_left < 0) || (ic_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		df_dcoord = meta::get_x_force(coord, d_coord, coord_min, potential_grid);
	}

    // df / dr = df / dcoord * dcoord / dr_ij * dr_ij / dr
    return (number) sign * r_vec * (df_dcoord * dcoord_dr / r);
}

number LTCoordination::potential(llint step, LR_vector &pos) {
    number coord = Utils::clamp(_coordination(), coord_min, coord_max); // ensure we are within the grid limits
    int ic_left = std::floor((coord - coord_min) / d_coord);
	int ic_right = ic_left + 1;

	number my_potential = 0;
	if((ic_left < 0) || (ic_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}
	else {
		my_potential = meta::interpolate_potential(coord, d_coord, coord_min, potential_grid);
	}

    return 10 * my_potential / all_pairs.size();
}

double LTCoordination::_coordination() {
    number coordination = 0.0;
    for(auto &pair : all_pairs) {
        number r = std::sqrt(CONFIG_INFO->box->sqr_min_image_distance(pair.first->pos, pair.second->pos));
        coordination += 1.0 / (1.0 + std::pow((r - d0) / r0, n));
    }

    return coordination;
}

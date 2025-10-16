/*
 * LTCoordination.cpp
 *
 * Created on: 10/16/2025
 *     Author: Lorenzo
*/

#include "LTCoordination.h"

#include "meta_utils.h"
#include "../../Utilities/OrderParameters.h"

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

    auto all_pairs = op.get_hb_particle_list();
    std::vector<int> p_indices;
    p_indices.reserve(all_pairs.size() * 2);
    for(auto &pair : all_pairs) {
        p_indices.push_back(pair.first);
        p_indices.push_back(pair.second);

        p1_ptr.push_back(CONFIG_INFO->particles()[pair.first]);
        p2_ptr.push_back(CONFIG_INFO->particles()[pair.second]);
    }

    // sorting is not necessary, but nice if we want to print all indices for debugging purposes
    // and also convenient to check for duplicates
    std::sort(p_indices.begin(), p_indices.end());

    // here we look for duplicates
    auto it = std::adjacent_find(p_indices.begin(), p_indices.end());
    if(it != p_indices.end()) {
        throw oxDNAException("LTCoordination: duplicate particle index found in the op_file: %d", *it);
    }

    // build the grid
    getInputNumber(&inp, "xmin", &xmin, 1);
	getInputNumber(&inp, "xmax", &xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &N_grid, 1);
	potential_grid.reserve(N_grid);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	potential_grid = meta::split_to_numbers(potential_string, ",");

	dX = (xmax - xmin) / (N_grid - 1.0);

    return std::make_tuple(p_indices, "LTCoordination force");
}

LR_vector LTCoordination::value(llint step, LR_vector &pos) {
    number coord = _coordination();
    int ix_left = std::floor((coord - xmin) / dX);
	int ix_right = ix_left + 1;

    
}

number LTCoordination::potential(llint step, LR_vector &pos) {
    number coord = _coordination();
    int ix_left = std::floor((coord - xmin) / dX);
	int ix_right = ix_left + 1;

	number my_potential = 0;
	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}
	else {
		my_potential = meta::interpolate_potential(coord, dX, xmin, potential_grid);
	}

    number total_factor = p1_ptr.size() + p2_ptr.size();
    return my_potential / total_factor;
}

double LTCoordination::_coordination() {
    return 0;
}

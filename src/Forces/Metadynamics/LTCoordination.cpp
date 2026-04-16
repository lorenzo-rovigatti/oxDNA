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
#include "../../Interactions/BaseInteraction.h"
#include "../../Interactions/DNA2Interaction.h"

#include <iostream>

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

    settings.get_settings(inp);

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

    std::string msg = Utils::sformat("LTCoordination force (d0 = %lf, r0 = %lf, n = %d, coord_min = %lf, coord_max = %lf, N_grid = %d)", settings.d0, settings.r0, settings.n, coord_min, coord_max, N_grid);
    return std::make_tuple(p_indices, msg);
}

LR_vector LTCoordination::force(llint step, LR_vector &pos) {
    auto pair = all_pairs[particle_to_pair_index[_current_particle]];
    int sign = (pair.first == _current_particle) ? -1 : 1;

    LR_vector r_vec = meta::distance(pair);
    number r = r_vec.module();

    number coord = Utils::clamp(meta::coordination(settings, all_pairs), coord_min, coord_max); // ensure we are within the grid limits
    int ic_left = std::floor((coord - coord_min) / d_coord);
	int ic_right = ic_left + 1;

    number df_dcoord = 0;
    // above we ensured coord is within limits, so this should never happen
	if((ic_left < 0) || (ic_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		df_dcoord = meta::get_x_force(coord, d_coord, coord_min, potential_grid);
	}

    // Numerically compute dcoord/dr for any coordination definition
    number dcoord_dr = _dcoord_dr(pair);

    LR_vector my_F = -sign * r_vec * (dcoord_dr / r);
    printf("FORCE pair: %d %d, idx: %d, %g (%g %g %g)\n", pair.first->index, pair.second->index, _current_particle->index, dcoord_dr, my_F.x, my_F.y, my_F.z);

    // df / dr = df / dcoord * dcoord / dr_ij * dr_ij / dr
    return (number) sign * r_vec * (df_dcoord * dcoord_dr / r);
}

LR_vector LTCoordination::torque(llint step, LR_vector &pos) {
    number coord = Utils::clamp(meta::coordination(settings, all_pairs), coord_min, coord_max);
    int ic_left = std::floor((coord - coord_min) / d_coord);

    number df_dcoord = 0;
    if((ic_left >= 0) && (ic_left + 1 <= N_grid - 1)) {
        df_dcoord = meta::get_x_force(coord, d_coord, coord_min, potential_grid);
    }

    // Part 1: Torque from off-center force application
    LR_vector F = force(step, pos);
    LR_vector torque = _current_particle->int_centers[DNANucleotide::BASE].cross(F);

    // Part 2: Direct torque from orientation-dependent coordination
    LR_vector dcoord_dtheta = _dcoord_dtheta();
    torque += dcoord_dtheta * df_dcoord;

    return torque;
}

number LTCoordination::potential(llint step, LR_vector &pos) {
    number coord = Utils::clamp(meta::coordination(settings, all_pairs), coord_min, coord_max); // ensure we are within the grid limits
    int ic_left = std::floor((coord - coord_min) / d_coord);
	int ic_right = ic_left + 1;

	number my_potential = 0;
	if((ic_left < 0) || (ic_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}
	else {
		my_potential = meta::interpolate_potential(coord, d_coord, coord_min, potential_grid);
	}

    return my_potential / (2.0 * all_pairs.size()); // divide by 2 to avoid double counting
}

number LTCoordination::_dcoord_dr(std::pair<BaseParticle*, BaseParticle*> &pair) {
    static number delta = 1e-6;
    
    LR_vector r_vec = meta::distance(pair);
    number r = r_vec.module();

    LR_vector dr_hat = r_vec / r;
    LR_vector old_pos = _current_particle->pos;
    // Determine the sign based on whether the current particle is the first or second in the pair
    int sign = (pair.first == _current_particle) ? 1 : -1;
    
    // Central difference - only compute contribution from this pair (particle belongs to single pair)
    _current_particle->pos = old_pos + dr_hat * (sign * delta);
    number contrib_plus = meta::get_pair_contribution(settings, pair);
    
    _current_particle->pos = old_pos - dr_hat * (sign * delta);
    number contrib_minus = meta::get_pair_contribution(settings, pair);
    
    _current_particle->pos = old_pos;  // Restore
    
    return (contrib_plus - contrib_minus) / (2.0 * delta);
}

LR_vector LTCoordination::_dcoord_dtheta() {
    static number delta = 1e-1;
    static number sin_delta = sin(delta);
    static number cos_delta = cos(delta);
    
    // Get the pair this particle belongs to
    auto pair = all_pairs[particle_to_pair_index[_current_particle]];
    
    LR_matrix old_orient = _current_particle->orientation;
    LR_matrix old_orientT = _current_particle->orientationT;
    
    LR_vector dcoord_dtheta(0, 0, 0);
    number contrib_plus, contrib_minus;
    
    // Central difference: Perturb rotation around x-axis (positive)
    _current_particle->orientation.v2 = old_orient.v2 * cos_delta + old_orient.v3 * sin_delta;
    _current_particle->orientation.v3 = -old_orient.v2 * sin_delta + old_orient.v3 * cos_delta;
    _current_particle->orientationT = _current_particle->orientation.get_transpose();
    _current_particle->set_positions();
    contrib_plus = meta::get_pair_contribution(settings, pair);
    
    // Central difference: Perturb rotation around x-axis (negative)
    _current_particle->orientation = old_orient;
    _current_particle->orientation.v2 = old_orient.v2 * cos_delta - old_orient.v3 * sin_delta;
    _current_particle->orientation.v3 = old_orient.v2 * sin_delta + old_orient.v3 * cos_delta;
    _current_particle->orientationT = _current_particle->orientation.get_transpose();
    _current_particle->set_positions();
    contrib_minus = meta::get_pair_contribution(settings, pair);
    dcoord_dtheta.x = (contrib_plus - contrib_minus) / (2.0 * delta);
    
    // Central difference: Perturb rotation around y-axis (positive)
    _current_particle->orientation = old_orient;
    _current_particle->orientation.v1 = old_orient.v1 * cos_delta - old_orient.v3 * sin_delta;
    _current_particle->orientation.v3 = old_orient.v1 * sin_delta + old_orient.v3 * cos_delta;
    _current_particle->orientationT = _current_particle->orientation.get_transpose();
    _current_particle->set_positions();
    contrib_plus = meta::get_pair_contribution(settings, pair);
    
    // Central difference: Perturb rotation around y-axis (negative)
    _current_particle->orientation = old_orient;
    _current_particle->orientation.v1 = old_orient.v1 * cos_delta + old_orient.v3 * sin_delta;
    _current_particle->orientation.v3 = -old_orient.v1 * sin_delta + old_orient.v3 * cos_delta;
    _current_particle->orientationT = _current_particle->orientation.get_transpose();
    _current_particle->set_positions();
    contrib_minus = meta::get_pair_contribution(settings, pair);
    dcoord_dtheta.y = (contrib_plus - contrib_minus) / (2.0 * delta);
    
    // Central difference: Perturb rotation around z-axis (positive)
    _current_particle->orientation = old_orient;
    _current_particle->orientation.v1 = old_orient.v1 * cos_delta + old_orient.v2 * sin_delta;
    _current_particle->orientation.v2 = -old_orient.v1 * sin_delta + old_orient.v2 * cos_delta;
    _current_particle->orientationT = _current_particle->orientation.get_transpose();
    _current_particle->set_positions();
    contrib_plus = meta::get_pair_contribution(settings, pair);
    
    // Central difference: Perturb rotation around z-axis (negative)
    _current_particle->orientation = old_orient;
    _current_particle->orientation.v1 = old_orient.v1 * cos_delta - old_orient.v2 * sin_delta;
    _current_particle->orientation.v2 = old_orient.v1 * sin_delta + old_orient.v2 * cos_delta;
    _current_particle->orientationT = _current_particle->orientation.get_transpose();
    _current_particle->set_positions();
    contrib_minus = meta::get_pair_contribution(settings, pair);
    dcoord_dtheta.z = (contrib_plus - contrib_minus) / (2.0 * delta);
    
    // Restore original orientation
    _current_particle->orientation = old_orient;
    _current_particle->orientationT = old_orientT;
    _current_particle->set_positions();

    printf("TORQUE pair: %d %d, idx: %d, %g %g %g\n", pair.first->index, pair.second->index, _current_particle->index, dcoord_dtheta.x, dcoord_dtheta.y, dcoord_dtheta.z);
    
    return dcoord_dtheta;
}

// LR_vector LTCoordination::_dcoord_dtheta() {
//     static number delta = 1e-1;
//     static number sin_delta = sin(delta);
//     static number cos_delta = cos(delta);
    
//     // Get the pair this particle belongs to
//     auto pair = all_pairs[particle_to_pair_index[_current_particle]];
//     number contrib_0 = meta::get_pair_contribution(settings, pair);
    
//     LR_matrix old_orient = _current_particle->orientation;
//     LR_matrix old_orientT = _current_particle->orientationT;
    
//     LR_vector dcoord_dtheta(0, 0, 0);
    
//     // Perturb rotation around x-axis
//     _current_particle->orientation.v2 = old_orient.v2 * cos_delta + old_orient.v3 * sin_delta;
//     _current_particle->orientation.v3 = -old_orient.v2 * sin_delta + old_orient.v3 * cos_delta;
//     _current_particle->orientationT = _current_particle->orientation.get_transpose();
//     _current_particle->set_positions();
//     dcoord_dtheta.x = (meta::get_pair_contribution(settings, pair) - contrib_0) / delta;
    
//     // Perturb rotation around y-axis
//     _current_particle->orientation = old_orient;
//     _current_particle->orientation.v1 = old_orient.v1 * cos_delta - old_orient.v3 * sin_delta;
//     _current_particle->orientation.v3 = old_orient.v1 * sin_delta + old_orient.v3 * cos_delta;
//     _current_particle->orientationT = _current_particle->orientation.get_transpose();
//     _current_particle->set_positions();
//     dcoord_dtheta.y = (meta::get_pair_contribution(settings, pair) - contrib_0) / delta;
    
//     // Perturb rotation around z-axis
//     _current_particle->orientation = old_orient;
//     _current_particle->orientation.v1 = old_orient.v1 * cos_delta + old_orient.v2 * sin_delta;
//     _current_particle->orientation.v2 = -old_orient.v1 * sin_delta + old_orient.v2 * cos_delta;
//     _current_particle->orientationT = _current_particle->orientation.get_transpose();
//     _current_particle->set_positions();
//     dcoord_dtheta.z = (meta::get_pair_contribution(settings, pair) - contrib_0) / delta;
    
//     // Restore original orientation
//     _current_particle->orientation = old_orient;
//     _current_particle->orientationT = old_orientT;
//     _current_particle->set_positions();

//     printf("TORQUE pair: %d %d, idx: %d, %g %g %g\n", pair.first->index, pair.second->index, _current_particle->index, dcoord_dtheta.x, dcoord_dtheta.y, dcoord_dtheta.z);
    
//     return dcoord_dtheta;
// }
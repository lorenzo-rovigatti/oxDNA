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
        throw oxDNAException("LTCoordination: duplicate particle index found in the op_file: %d. LTCoordination assumes each particle appears only once", *it);
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

    // initialise rotation matrices for torque calculation
    LR_vector axes[3] = {LR_vector(1, 0, 0), LR_vector(0, 1, 0), LR_vector(0, 0, 1)};
    for(int i = 0; i < 3; i++) {
        _rot_matrices[i][0] = Utils::get_rotation_matrix_from_axis_angle(axes[i], delta_T);
        _rot_matrices[i][1] = Utils::get_rotation_matrix_from_axis_angle(axes[i], -delta_T);;
    }

    std::string msg = Utils::sformat("LTCoordination force (d0 = %lf, r0 = %lf, n = %d, coord_min = %lf, coord_max = %lf, N_grid = %d)", settings.d0, settings.r0, settings.n, coord_min, coord_max, N_grid);
    return std::make_tuple(p_indices, msg);
}

LR_vector LTCoordination::force(llint step, LR_vector &pos) {
    number coord = _coordination(step);
    int ic_left = std::floor((coord - coord_min) / d_coord);
	int ic_right = ic_left + 1;

    number df_dcoord = 0;
    // above we ensured coord is within limits, so this should never happen
	if((ic_left < 0) || (ic_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
        // note that this is not striclty only df/dcoord, since it also has a minux sign (hence the "get_x_force" name)
		df_dcoord = meta::get_x_force(coord, d_coord, coord_min, potential_grid);
	}

    LR_vector dcoord_dpos = _dcoord_dpos();

    // DNANucleotide p = *dynamic_cast<DNANucleotide *>(_current_particle == pair.first ? pair.first : pair.second);
    // DNANucleotide q = *dynamic_cast<DNANucleotide *>(_current_particle == pair.first ? pair.second : pair.first);
    // p.force = LR_vector();
    // q.force = LR_vector();
    // CONFIG_INFO->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, &p, &q, true, true);
    // printf("FORCE %d\n\t%g %g %g\n", p.index, p.force.x, p.force.y, p.force.z);
    // printf("\t%g %g %g\n", dcoord_dpos.x, dcoord_dpos.y, dcoord_dpos.z);

    return dcoord_dpos * df_dcoord;
}

LR_vector LTCoordination::torque(llint step, LR_vector &pos) {
    number coord = _coordination(step);
    int ic_left = std::floor((coord - coord_min) / d_coord);

    number df_dcoord = 0;
    if((ic_left >= 0) && (ic_left + 1 <= N_grid - 1)) {
        // note that this is not striclty only df/dcoord, since it also has a minux sign (hence the "get_x_force" name)
        df_dcoord = meta::get_x_force(coord, d_coord, coord_min, potential_grid);
    }

    LR_vector dcoord_dtheta = _dcoord_dtheta();

    // auto pair = all_pairs[particle_to_pair_index[_current_particle]];
    // DNANucleotide p = *dynamic_cast<DNANucleotide *>(_current_particle == pair.first ? pair.first : pair.second);
    // DNANucleotide q = *dynamic_cast<DNANucleotide *>(_current_particle == pair.first ? pair.second : pair.first);
    // p.force = LR_vector();
    // q.force = LR_vector();
    // LR_vector pref_torque = _current_particle->orientationT * dcoord_dtheta;
    // CONFIG_INFO->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, &p, &q, true, true);
    // printf("TORQUE %d\n\t%g %g %g\n", p.index, p.torque.x, p.torque.y, p.torque.z);
    // printf("\t%g %g %g\n", pref_torque.x, pref_torque.y, pref_torque.z);

    return dcoord_dtheta * df_dcoord;
}

number LTCoordination::potential(llint step, LR_vector &pos) {
    number coord = _coordination(step);
    int ic_left = std::floor((coord - coord_min) / d_coord);
	int ic_right = ic_left + 1;

	number my_potential = 0;
	if((ic_left < 0) || (ic_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}
	else {
		my_potential = meta::interpolate_potential(coord, d_coord, coord_min, potential_grid);
        // printf("coord = %lf, index = %d, potential = %lf\n", coord, ic_left, my_potential);
	}

    return my_potential / (2.0 * all_pairs.size()); // divide by 2 to avoid double counting
}

LR_vector LTCoordination::_dcoord_dpos() {
    auto pair = all_pairs[particle_to_pair_index[_current_particle]];
    LR_vector old_pos = _current_particle->pos;
    LR_vector grad;

    for(int i = 0; i < 3; i++) {
        LR_vector shift;
        shift[i] = delta_F;

        _current_particle->pos = old_pos + shift;
        number plus = meta::get_pair_contribution(settings, pair);

        _current_particle->pos = old_pos - shift;
        number minus = meta::get_pair_contribution(settings, pair);

        grad[i] = (plus - minus) / (2.0 * delta_F);
    }
    _current_particle->pos = old_pos;
    return grad;
}

LR_vector LTCoordination::_dcoord_dtheta() {
    // Get the pair this particle belongs to
    auto pair = all_pairs[particle_to_pair_index[_current_particle]];
    
    LR_matrix old_orient = _current_particle->orientation;
    LR_matrix old_orientT = _current_particle->orientationT;
    
    LR_vector dcoord_dtheta(0, 0, 0);
    number contrib_plus, contrib_minus;

    // here we use a central difference scheme, perturbing the orientation around each axis 
    // *in the lab reference frame* in both directions and computing the contribution to the 
    // coordination from the pair this particle belongs to
    for(int i = 0; i < 3; i++) {
        // positive rotation
        _current_particle->orientation = _rot_matrices[i][0] * old_orient;
        _current_particle->orientationT = _current_particle->orientation.get_transpose();
        _current_particle->set_positions();
        contrib_plus = meta::get_pair_contribution(settings, pair);

        // negative rotation
        _current_particle->orientation = _rot_matrices[i][1] * old_orient;
        _current_particle->orientationT = _current_particle->orientation.get_transpose();
        _current_particle->set_positions();
        contrib_minus = meta::get_pair_contribution(settings, pair);

        dcoord_dtheta[i] = (contrib_plus - contrib_minus) / (2.0 * delta_T);
    }
    
    // Restore original orientation
    _current_particle->orientation = old_orient;
    _current_particle->orientationT = old_orientT;
    _current_particle->set_positions();
    
    return dcoord_dtheta;
}

number LTCoordination::_coordination(llint step) {
    if(step != _last_step_calculated) {
        _current_coordination = Utils::clamp(meta::coordination(settings, all_pairs), coord_min, coord_max);
        _last_step_calculated = step;
    }
    return _current_coordination;
}

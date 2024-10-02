#include "LTCOMTrap.h"

#include "meta_utils.h"
#include "../../Particles/BaseParticle.h"

#include <string>
#include <iostream>

inline number interpolatePotential(const number my_x, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	number x_left = dX * std::floor(my_x / dX);
	number x_right = x_left + dX;
	number ix_left = std::floor((my_x - xmin) / dX);
	number ix_right = ix_left + 1;
	number f11 = potential_grid[ix_left];
	number f21 = potential_grid[ix_right];
	number fx = (x_right - my_x) / dX * f11 + (my_x - x_left) / dX * f21;
	return fx;
}

inline number get_x_force(const number x, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	number ix_left = std::floor((x - xmin) / dX);
	number ix_right = ix_left + 1;
	return -(potential_grid[ix_right] - potential_grid[ix_left]) / dX;
}

LTCOMTrap::LTCOMTrap() :
				BaseForce() {
	xmin = 0;
	xmax = 10;
	N_grid = 0;
	dX = 0.1;
}

std::tuple<std::vector<int>, std::string> LTCOMTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::tie(_p1a, _p1a_ptr) = meta::get_particle_lists(inp, "p1a", CONFIG_INFO->particles(), "LTCOMTrap p1a");
	std::tie(_p2a, _p2a_ptr) = meta::get_particle_lists(inp, "p2a", CONFIG_INFO->particles(), "LTCOMTrap p2a");

	getInputBool(&inp, "PBC", &PBC, 0);
	getInputInt(&inp, "mode", &_mode, 1);

	if(_mode != 1 and _mode != 2) {
		throw oxDNAException("LTCOMTrap: unsupported mode '%d' (should be '1' or '2')", _mode);
	}

	getInputNumber(&inp, "xmin", &xmin, 1);
	getInputNumber(&inp, "xmax", &xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &N_grid, 1);
	potential_grid.reserve(N_grid);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	potential_grid = meta::split_to_numbers(potential_string, ",");

	dX = (xmax - xmin) / (N_grid - 1.0);

	std::string description = Utils::sformat("LTCOMTrap force with mode = %d", _mode);
	if(_mode == 1) {
		return std::make_tuple(_p1a, description);
	}
	else {
		return std::make_tuple(_p2a, description);
	}
}

LR_vector LTCOMTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC)
		return CONFIG_INFO->box->min_image(u, v);
	else
		return v - u;
}

LR_vector LTCOMTrap::value(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);

	LR_vector dra = _distance(p2a_vec, p1a_vec);

	number my_x = dra.module();

	int ix_left = std::floor((my_x - xmin) / dX);
	int ix_right = ix_left + 1;

	number meta_Fx = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		meta_Fx = get_x_force(my_x, dX, xmin, potential_grid);
	}

	const LR_vector accumulated_force = dra * (meta_Fx / dra.module());

	if(_mode == 1) {
		return accumulated_force / (number) _p1a_ptr.size();
	}
	else {
		return accumulated_force / (-1 * (number) _p2a_ptr.size());
	}
}

number LTCOMTrap::potential(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);

	LR_vector dra = _distance(p2a_vec, p1a_vec);

	number x = dra.module();

	int ix_left = std::floor((x - xmin) / dX);
	int ix_right = ix_left + 1;

	number my_potential = 0;
	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		my_potential = interpolatePotential(x, dX, xmin, potential_grid);
	}

	if(_mode == 1) {
		return my_potential / (number) _p1a_ptr.size();
	}
	else {
		return my_potential / ((number) _p2a_ptr.size());
	}

	number total_factor = _p1a_ptr.size() + _p2a_ptr.size();

	return my_potential / total_factor;
}

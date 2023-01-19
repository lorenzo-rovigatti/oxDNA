#include "LTCOMAngleTrap.h"

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

LTCOMAngleTrap::LTCOMAngleTrap() :
				BaseForce() {
	xmin = 0;
	xmax = 10;
	N_grid = 0;
	dX = 0.1;
}

std::tuple<std::vector<int>, std::string> LTCOMAngleTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::tie(_p1a, _p1a_ptr) = meta::get_particle_lists(inp, "p1a", CONFIG_INFO->particles(), "GaussTrap p1a");
	std::tie(_p2a, _p2a_ptr) = meta::get_particle_lists(inp, "p2a", CONFIG_INFO->particles(), "GaussTrap p2a");
	std::tie(_p3a, _p3a_ptr) = meta::get_particle_lists(inp, "p3a", CONFIG_INFO->particles(), "GaussTrap p3a");

	getInputInt(&inp, "mode", &_mode, 1);
	if(_mode != 1 and _mode != 2 and _mode != 3) {
		throw oxDNAException("LTCOMAngleTrap: unsupported mode '%d' (should be '1', '2' or '3')", _mode);
	}

	getInputBool(&inp, "PBC", &PBC, 0);

	getInputNumber(&inp, "xmin", &xmin, 1);
	getInputNumber(&inp, "xmax", &xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &N_grid, 1);
	potential_grid.reserve(N_grid);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	potential_grid = meta::split_to_numbers(potential_string, ",");

	dX = (xmax - xmin) / (N_grid - 1.0);

	std::string description = Utils::sformat("LTCOMAngleTrap force with mode = %d", _mode);
	if(_mode == 1) {
		return std::make_tuple(_p1a, description);
	}
	else if(_mode == 2) {
		return std::make_tuple(_p2a, description);
	}
	else {
		return std::make_tuple(_p3a, description);
	}
}

LR_vector LTCOMAngleTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC)
		return CONFIG_INFO->box->min_image(u, v);
	else
		return v - u;
}

LR_vector LTCOMAngleTrap::value(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);
	LR_vector p3a_vec = meta::particle_list_com(_p3a_ptr, CONFIG_INFO->box);

	LR_vector dra1 = _distance(p2a_vec, p1a_vec);
	LR_vector dra2 = _distance(p2a_vec, p3a_vec);

	LR_vector dra1_normed = dra1 / dra1.module();
	LR_vector dra2_normed = dra2 / dra2.module();

	number dot_product = dra1_normed * dra2_normed;
	number angle = std::acos(dot_product);
	int ix_left = std::floor((angle - xmin) / dX);
	int ix_right = ix_left + 1;

	number xforce = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		xforce = get_x_force(angle, dX, xmin, potential_grid);
	}

	number prefactor = -xforce / std::pow(1 - dot_product, 0.5);
	if(_mode == 1) {
		number r1_factor = -dot_product / dra1.module();
		number r2_factor = 1 / dra1.module();
		return ((dra1_normed * r1_factor) + (dra2_normed * r2_factor)) * prefactor / (number) _p1a_ptr.size();
	}
	if(_mode == 2) {
		number r1_factor_A = -dot_product / dra1.module();
		number r2_factor_A = 1 / dra1.module();
		number r2_factor_B = -dot_product / dra2.module();
		number r1_factor_B = 1 / dra2.module();
		return ((dra1_normed * (-r1_factor_A - r1_factor_B)) + (dra2_normed * (-r2_factor_A - r2_factor_B))) * prefactor / (number) _p2a_ptr.size();
	}
	else {
		number r2_factor = -dot_product / dra2.module();
		number r1_factor = 1 / dra2.module();
		return ((dra1_normed * r1_factor) + (dra2_normed * r2_factor)) * prefactor / (number) _p3a_ptr.size();
	}
}

number LTCOMAngleTrap::potential(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);
	LR_vector p3a_vec = meta::particle_list_com(_p3a_ptr, CONFIG_INFO->box);

	LR_vector dra1 = _distance(p2a_vec, p1a_vec);
	LR_vector dra2 = _distance(p2a_vec, p3a_vec);

	LR_vector dra1_normed = dra1 / dra1.module();
	LR_vector dra2_normed = dra2 / dra2.module();

	number dot_product = dra1_normed * dra2_normed;
	number angle = std::acos(dot_product);
	int ix_left = std::floor((angle - xmin) / dX);
	int ix_right = ix_left + 1;

	number my_potential = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		my_potential = interpolatePotential(angle, dX, xmin, potential_grid);
	}

	if(_mode == 1) {
		return my_potential / (number) _p1a_ptr.size();
	}
	if(_mode == 2) {
		return my_potential / (number) _p2a_ptr.size();
	}
	else {
		return my_potential / (number) _p3a_ptr.size();
	}
}

#include "LTAtanCOMTrap.h"

#include "meta_utils.h"
#include "../../Particles/BaseParticle.h"

#include <tuple>
#include <string>
#include <iostream>

inline number interpolatePotential2(number a, number dX, number xmin, const std::vector<number> &potential_grid) {
	number a_left = dX * std::floor(a / dX);
	number a_right = a_left + dX;

	number ia_left = std::floor((a - xmin) / dX);
	number ia_right = ia_left + 1;

	number f1 = potential_grid[ia_left];
	number f2 = potential_grid[ia_right];

	number f = (a_right - a) / dX * f1 + (a - a_left) / dX * f2;
	return f;
}

inline number get_x_force2(const number x, const number y, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	const number delta = dX / 5.;
	const number ratio = y / x;
	const number a = std::atan(ratio);
	const number Uprime = (interpolatePotential2(a + delta, dX, xmin, potential_grid) - interpolatePotential2(a - delta, dX, xmin, potential_grid)) / (2 * delta);
	const number prefactor = -(1. / (1 + std::pow(ratio, 2))) * (y / std::pow(x, 2)); // negative sign comes from derivative wrt denominator
	const number force = -prefactor * Uprime;
	return force;
}
inline number get_y_force2(const number x, const number y, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	const number delta = dX / 5.;
	const number ratio = y / x;
	const number a = std::atan(ratio);
	const number Uprime = (interpolatePotential2(a + delta, dX, xmin, potential_grid) - interpolatePotential2(a - delta, dX, xmin, potential_grid)) / (2 * delta);
	const number prefactor = (1. / (1 + std::pow(ratio, 2))) * (1 / x);
	const number force = -prefactor * Uprime;
	return force;
}

LTAtanCOMTrap::LTAtanCOMTrap() :
				BaseForce() {
	xmin = 0;
	xmax = 10;
	N_grid = 0;
	dX = 0.1;
}

std::tuple<std::vector<int>, std::string> LTAtanCOMTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::tie(_p1a, _p1a_ptr) = meta::get_particle_lists(inp, "p1a", CONFIG_INFO->particles(), "GaussTrapMetaRatio1D p1a");
	std::tie(_p2a, _p2a_ptr) = meta::get_particle_lists(inp, "p2a", CONFIG_INFO->particles(), "GaussTrapMetaRatio1D p2a");
	std::tie(_p1b, _p1b_ptr) = meta::get_particle_lists(inp, "p1b", CONFIG_INFO->particles(), "GaussTrapMetaRatio1D p1b");
	std::tie(_p2b, _p2b_ptr) = meta::get_particle_lists(inp, "p2b", CONFIG_INFO->particles(), "GaussTrapMetaRatio1D p2b");

	getInputInt(&inp, "mode", &_mode, 1);
	if(_mode != 1 and _mode != 2 and _mode != 3 and _mode != 4) {
		throw oxDNAException("LTAtanCOMTrap: unsupported mode '%d' (should be '1', '2', '3' or '4')", _mode);
	}

	getInputBool(&inp, "PBC", &PBC, 0);

	getInputNumber(&inp, "xmin", &xmin, 1);
	getInputNumber(&inp, "xmax", &xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &N_grid, 1);
	potential_grid.reserve(N_grid);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	potential_grid = meta::split_to_numbers(potential_string, ",");

	dX = (xmax - xmin) / (N_grid - 1);

	std::string description = Utils::sformat("LTAtanCOMTrap force with mode = %d", _mode);
	if(_mode == 1) {
		return std::make_tuple(_p1a, description);
	}
	else if(_mode == 2) {
		return std::make_tuple(_p2a, description);
	}
	else if(_mode == 3) {
		return std::make_tuple(_p1b, description);
	}
	else {
		return std::make_tuple(_p2b, description);
	}
	std::tuple<std::vector<int>, std::string>();
}

LR_vector LTAtanCOMTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC)
		return CONFIG_INFO->box->min_image(u, v);
	else
		return v - u;
}

LR_vector LTAtanCOMTrap::value(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);
	LR_vector p1b_vec = meta::particle_list_com(_p1b_ptr, CONFIG_INFO->box);
	LR_vector p2b_vec = meta::particle_list_com(_p2b_ptr, CONFIG_INFO->box);

	LR_vector dra = _distance(p2a_vec, p1a_vec);
	LR_vector drb = _distance(p2b_vec, p1b_vec);

	number my_x = dra.module();
	number my_y = drb.module();

	number angle = std::atan(my_y / my_x);
	int ix_left = std::floor((angle - xmin) / dX);
	int ix_right = ix_left + 1;

	number meta_Fx = 0;
	number meta_Fy = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		meta_Fx = get_x_force2(my_x, my_y, dX, xmin, potential_grid);
		meta_Fy = get_y_force2(my_x, my_y, dX, xmin, potential_grid);
		//std::cout << meta_Fx << std::endl;
		//std::cout << meta_Fy << std::endl;
	}

	const LR_vector accumulated_force_x = dra * (meta_Fx / dra.module());
	const LR_vector accumulated_force_y = drb * (meta_Fy / drb.module());

	if(_mode == 1) {
		return accumulated_force_x / _p1a.size();
	}
	else if(_mode == 2) {
		return (accumulated_force_x * -1) / _p2a.size();
	}
	else if(_mode == 3) {
		return accumulated_force_y / _p1b.size();
	}
	else {
		return (accumulated_force_y * -1) / _p2b.size();
	}
}

number LTAtanCOMTrap::potential(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);
	LR_vector p1b_vec = meta::particle_list_com(_p1b_ptr, CONFIG_INFO->box);
	LR_vector p2b_vec = meta::particle_list_com(_p2b_ptr, CONFIG_INFO->box);

	LR_vector dra = _distance(p2a_vec, p1a_vec);
	LR_vector drb = _distance(p2b_vec, p1b_vec);

	number my_x = dra.module();
	number my_y = drb.module();
	number angle = my_y / my_x;
	int ix_left = std::floor((angle - xmin) / dX);
	int ix_right = ix_left + 1;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	number my_potential = 0;

	const number ratio = my_y / my_x;
	const number a = std::atan(ratio);

	my_potential = interpolatePotential2(a, dX, xmin, potential_grid);

	if(_mode == 1) {
		return my_potential / _p1a.size();
	}
	else if(_mode == 2) {
		return my_potential / _p2a.size();
	}
	else if(_mode == 3) {
		return my_potential / _p1b.size();
	}
	else {
		return my_potential / _p2b.size();
	}
}

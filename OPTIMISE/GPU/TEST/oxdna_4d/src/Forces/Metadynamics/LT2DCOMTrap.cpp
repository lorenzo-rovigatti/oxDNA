#include "LT2DCOMTrap.h"

#include "meta_utils.h"
#include "../../Particles/BaseParticle.h"

#include <string>
#include <iostream>

number interpolatePotential2(number x, number y, number dX, number xmin, const std::vector<std::vector<number> > &potential_grid) {
	number x_left = dX * std::floor(x / dX);
	number y_left = dX * std::floor(y / dX);
	number x_right = x_left + dX;
	number y_right = y_left + dX;

	number ix_left = std::floor((x - xmin) / dX);
	number iy_left = std::floor((y - xmin) / dX);
	number ix_right = ix_left + 1;
	number iy_right = iy_left + 1;

	number f11 = potential_grid[ix_left][iy_left];
	number f12 = potential_grid[ix_left][iy_right];
	number f21 = potential_grid[ix_right][iy_left];
	number f22 = potential_grid[ix_right][iy_right];

	number fxy1 = (x_right - x) / dX * f11 + (x - x_left) / dX * f21;
	number fxy2 = (x_right - x) / dX * f12 + (x - x_left) / dX * f22;
	number local_potential = (y_right - y) / dX * fxy1 + (y - y_left) / dX * fxy2;
	return local_potential;
}

inline number get_x_force2(const number x, const number y, const number dX, const number xmin, const std::vector<std::vector<number> > &potential_grid) {
	const number delta = dX / 5.;
	return -(interpolatePotential2(x + delta, y, dX, xmin, potential_grid) - interpolatePotential2(x - delta, y, dX, xmin, potential_grid)) / (2 * delta);
}
inline number get_y_force2(const number x, const number y, const number dX, const number xmin, const std::vector<std::vector<number> > &potential_grid) {
	const number delta = dX / 5.;
	return -(interpolatePotential2(x, y + delta, dX, xmin, potential_grid) - interpolatePotential2(x, y - delta, dX, xmin, potential_grid)) / (2 * delta);
}

LT2DCOMTrap::LT2DCOMTrap() :
				BaseForce() {
	xmin = 0;
	xmax = 10;
	N_grid = 0;
	dX = 0.1;
}

std::tuple<std::vector<int>, std::string> LT2DCOMTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::tie(_p1a, _p1a_ptr) = meta::get_particle_lists(inp, "p1a", CONFIG_INFO->particles(), "GaussTrapMeta p1a");
	std::tie(_p2a, _p2a_ptr) = meta::get_particle_lists(inp, "p2a", CONFIG_INFO->particles(), "GaussTrapMeta p2a");
	std::tie(_p1b, _p1b_ptr) = meta::get_particle_lists(inp, "p1b", CONFIG_INFO->particles(), "GaussTrapMeta p1b");
	std::tie(_p2b, _p2b_ptr) = meta::get_particle_lists(inp, "p2b", CONFIG_INFO->particles(), "GaussTrapMeta p2b");

	getInputInt(&inp, "mode", &_mode, 1);
	if(_mode != 1 and _mode != 2 and _mode != 3 and _mode != 4) {
		throw oxDNAException("LT2DCOMTrap: unsupported mode '%d' (should be '1', '2', '3' or '4')", _mode);
	}

	getInputBool(&inp, "PBC", &PBC, 0);

	getInputNumber(&inp, "xmin", &xmin, 1);
	getInputNumber(&inp, "xmax", &xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &N_grid, 1);
	potential_grid.reserve(N_grid);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	for(auto &token : Utils::split(potential_string, '|')) {
		potential_grid.push_back(meta::split_to_numbers(token, ","));
	}

	dX = (xmax - xmin) / (N_grid - 1);

	std::string description = Utils::sformat("LT2DCOMTrap force with mode = %d", _mode);
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
}

LR_vector LT2DCOMTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC)
		return CONFIG_INFO->box->min_image(u, v);
	else
		return v - u;
}

LR_vector LT2DCOMTrap::value(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);
	LR_vector p1b_vec = meta::particle_list_com(_p1b_ptr, CONFIG_INFO->box);
	LR_vector p2b_vec = meta::particle_list_com(_p2b_ptr, CONFIG_INFO->box);

	LR_vector dra = _distance(p2a_vec, p1a_vec);
	LR_vector drb = _distance(p2b_vec, p1b_vec);

	number my_x = dra.module();
	number my_y = drb.module();

	int ix_left = std::floor((my_x - xmin) / dX);
	int ix_right = ix_left + 1;
	int iy_left = std::floor((my_y - xmin) / dX);
	int iy_right = iy_left + 1;

	number meta_Fx = 0;
	number meta_Fy = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}
	if((iy_left < 0) || (iy_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		meta_Fx = get_x_force2(my_x, my_y, dX, xmin, potential_grid);
		meta_Fy = get_y_force2(my_x, my_y, dX, xmin, potential_grid);
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

number LT2DCOMTrap::potential(llint step, LR_vector &pos) {
	LR_vector p1a_vec = meta::particle_list_com(_p1a_ptr, CONFIG_INFO->box);
	LR_vector p2a_vec = meta::particle_list_com(_p2a_ptr, CONFIG_INFO->box);
	LR_vector p1b_vec = meta::particle_list_com(_p1b_ptr, CONFIG_INFO->box);
	LR_vector p2b_vec = meta::particle_list_com(_p2b_ptr, CONFIG_INFO->box);

	LR_vector dra = _distance(p2a_vec, p1a_vec);
	LR_vector drb = _distance(p2b_vec, p1b_vec);

	number my_x = dra.module();
	number my_y = drb.module();

	int ix_left = std::floor((my_x - xmin) / dX);
	int ix_right = ix_left + 1;
	int iy_left = std::floor((my_y - xmin) / dX);
	int iy_right = iy_left + 1;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}
	if((iy_left < 0) || (iy_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	number my_potential = 0;

	my_potential = interpolatePotential2(my_x, my_y, dX, xmin, potential_grid);

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

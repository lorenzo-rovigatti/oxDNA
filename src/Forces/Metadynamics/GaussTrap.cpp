#include "GaussTrap.h"

#include "../../Particles/BaseParticle.h"
#include "../../Boxes/BaseBox.h"

#include <string>
#include <sstream>
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

GaussTrap::GaussTrap() :
				BaseForce() {
	xmin = 0;
	xmax = 10;
	N_grid = 0;
	dX = 0.1;
}

std::tuple<std::vector<int>, std::string> GaussTrap::init(input_file &inp, BaseBox *box_ptr) {
	std::string p1a_string;
	std::string p2a_string;

	getInputString(&inp, "p1a", p1a_string, 1);
	getInputString(&inp, "p2a", p2a_string, 1);
	getInputBool(&inp, "PBC", &PBC, 0);
	getInputInt(&inp, "mode", &_mode, 1);

	if(_mode != 1 and _mode != 2) {
		throw oxDNAException("GaussTrap: unsupported mode '%d' (should be '1' or '2')", _mode);
	}

	_p1a = Utils::get_particles_from_string(CONFIG_INFO->particles(), p1a_string, "GaussTrap p1a");
	for(auto p_idx : _p1a) {
		_p1a_ptr.push_back(CONFIG_INFO->particles()[p_idx]);
	}

	_p2a = Utils::get_particles_from_string(CONFIG_INFO->particles(), p2a_string, "GaussTrap p2a");
	for(auto p_idx : _p2a) {
		_p2a_ptr.push_back(CONFIG_INFO->particles()[p_idx]);
	}

	getInputNumber(&inp, "xmin", &xmin, 1);
	getInputNumber(&inp, "xmax", &xmax, 1); // these are both inclusive
	getInputInt(&inp, "N_grid", &N_grid, 1);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	for(auto token : Utils::split(potential_string, ',')) {
		potential_grid.push_back(std::stod(token));
	}

	_box_ptr = box_ptr;
	dX = (xmax - xmin) / (N_grid - 1.0);

	std::string description = Utils::sformat("GaussTrap force with mode = %d", _mode);
	if(_mode == 1) {
		return std::make_tuple(_p1a, description);
	}
	else {
		return std::make_tuple(_p2a, description);
	}
}

LR_vector GaussTrap::_distance(LR_vector u, LR_vector v) {
	if(this->PBC)
		return this->_box_ptr->min_image(u, v);
	else
		return v - u;
}

LR_vector GaussTrap::value(llint step, LR_vector &pos) {
	LR_vector p1a_vec;
	LR_vector p2a_vec;

	for(auto p : _p1a_ptr) {
		p1a_vec += _box_ptr->get_abs_pos(p);
	}
	p1a_vec = p1a_vec / (number) _p1a_ptr.size();

	for(auto p : _p2a_ptr) {
		p2a_vec += _box_ptr->get_abs_pos(p);
	}
	p2a_vec = p2a_vec / (number) _p2a_ptr.size();

	LR_vector dra = _distance(p2a_vec, p1a_vec);

	number my_x = dra.module();

	int ix_left = std::floor((my_x - this->xmin) / this->dX);
	int ix_right = ix_left + 1;

	number meta_Fx = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		meta_Fx = get_x_force(my_x, dX, xmin, potential_grid);
	}

	const LR_vector accumulated_force = dra * (meta_Fx / dra.module());

	if(this->_mode == 1) {
		return accumulated_force / (number) _p1a_ptr.size();
	}
	else {
		return accumulated_force / (-1 * (number) _p2a_ptr.size());
	}
}

number GaussTrap::potential(llint step, LR_vector &pos) {
	LR_vector p1a_vec = { 0, 0, 0 };
	LR_vector p2a_vec = { 0, 0, 0 };

	for(auto p : _p1a_ptr) {
		p1a_vec += _box_ptr->get_abs_pos(p);
	}
	p1a_vec = p1a_vec / (number) _p1a_ptr.size();

	for(auto p : _p2a_ptr) {
		p2a_vec += _box_ptr->get_abs_pos(p);
	}
	p2a_vec = p2a_vec / (number) _p2a_ptr.size();

	LR_vector dra = this->_distance(p2a_vec, p1a_vec);

	number x = dra.module();

	int ix_left = std::floor((x - this->xmin) / this->dX);
	int ix_right = ix_left + 1;

	number my_potential = 0;
	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		my_potential = interpolatePotential(x, dX, xmin, potential_grid);
	}

	number total_factor = _p1a_ptr.size() + _p2a_ptr.size();

	return my_potential / total_factor;
}

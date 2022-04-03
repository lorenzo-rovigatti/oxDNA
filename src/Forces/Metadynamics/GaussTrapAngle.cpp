#include "GaussTrapAngle.h"
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

GaussTrapAngle::GaussTrapAngle() :
				BaseForce() {
	xmin = 0;
	xmax = 10;
	N_grid = 0;
	dX = 0.1;
}

std::tuple<std::vector<int>, std::string> GaussTrapAngle::init(input_file &inp, BaseBox *box_ptr) {
	std::string p1a_string;
	std::string p2a_string;
	std::string p3a_string;

	getInputString(&inp, "p1a", p1a_string, 1);
	getInputString(&inp, "p2a", p2a_string, 1);
	getInputString(&inp, "p3a", p3a_string, 1);

	_p1a = Utils::get_particles_from_string(CONFIG_INFO->particles(), p1a_string, "GaussTrap p1a");
	for(auto p_idx : _p1a) {
		_p1a_ptr.push_back(CONFIG_INFO->particles()[p_idx]);
	}

	_p2a = Utils::get_particles_from_string(CONFIG_INFO->particles(), p2a_string, "GaussTrap p2a");
	for(auto p_idx : _p2a) {
		_p2a_ptr.push_back(CONFIG_INFO->particles()[p_idx]);
	}

	_p3a = Utils::get_particles_from_string(CONFIG_INFO->particles(), p2a_string, "GaussTrap p3a");
	for(auto p_idx : _p3a) {
		_p3a_ptr.push_back(CONFIG_INFO->particles()[p_idx]);
	}

	getInputInt(&inp, "mode", &_mode, 1);
	if(_mode != 1 and _mode != 2 and _mode != 3) {
		throw oxDNAException("GaussTrapAngle: unsupported mode '%d' (should be '1', '2' or '3')", _mode);
	}

	getInputBool(&inp, "PBC", &PBC, 0);

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

	std::string description = Utils::sformat("GaussTrapAngle force with mode = %d", _mode);
	if(this->_mode == 1) {
		return std::make_tuple(_p1a, description);
	}
	else if(this->_mode == 2) {
		return std::make_tuple(_p2a, description);
	}
	else {
		return std::make_tuple(_p3a, description);
	}
}

LR_vector GaussTrapAngle::_distance(LR_vector u, LR_vector v) {
	if(this->PBC)
		return this->_box_ptr->min_image(u, v);
	else
		return v - u;
}

LR_vector GaussTrapAngle::value(llint step, LR_vector &pos) {
	LR_vector p1a_vec = { 0, 0, 0 };
	LR_vector p2a_vec = { 0, 0, 0 };
	LR_vector p3a_vec = { 0, 0, 0 };

	for(auto p : _p1a_ptr) {
		p1a_vec += this->_box_ptr->get_abs_pos(p);
	}
	p1a_vec = p1a_vec / (number) _p1a_ptr.size();

	for(auto p : _p2a_ptr) {
		p2a_vec += this->_box_ptr->get_abs_pos(p);
	}
	p2a_vec = p2a_vec / (number) _p2a_ptr.size();

	for(auto p : _p3a_ptr) {
		p3a_vec += this->_box_ptr->get_abs_pos(p);
	}
	p3a_vec = p3a_vec / (number) _p3a_ptr.size();

	LR_vector dra1 = this->_distance(p2a_vec, p1a_vec);
	LR_vector dra2 = this->_distance(p2a_vec, p3a_vec);

	LR_vector dra1_normed = dra1 / dra1.module();
	LR_vector dra2_normed = dra2 / dra2.module();

	number dot_product = dra1_normed * dra2_normed;
	number angle = std::acos(dot_product);
	int ix_left = std::floor((angle - this->xmin) / this->dX);
	int ix_right = ix_left + 1;

	number xforce = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		xforce = get_x_force(angle, dX, xmin, potential_grid);
	}

	// why is this necessary?
	number prefactor = -xforce / std::pow(1 - dot_product, 0.5);

	// this isn't fucking working
	if(this->_mode == 1) {
		number r1_factor = -dot_product / dra1.module();
		number r2_factor = 1 / dra1.module();
		return ((dra1_normed * r1_factor) + (dra2_normed * r2_factor)) * prefactor / (number) _p1a_ptr.size();
	}
	if(this->_mode == 2) {
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

number GaussTrapAngle::potential(llint step, LR_vector &pos) {
	LR_vector p1a_vec = { 0, 0, 0 };
	LR_vector p2a_vec = { 0, 0, 0 };
	LR_vector p3a_vec = { 0, 0, 0 };

	for(auto p : _p1a_ptr) {
		p1a_vec += _box_ptr->get_abs_pos(p);
	}
	p1a_vec = p1a_vec / (number) _p1a_ptr.size();

	for(auto p : _p2a_ptr) {
		p2a_vec += _box_ptr->get_abs_pos(p);
	}
	p2a_vec = p2a_vec / (number) _p2a_ptr.size();

	for(auto p : _p3a_ptr) {
		p3a_vec += _box_ptr->get_abs_pos(p);
	}
	p3a_vec = p3a_vec / (number) _p3a_ptr.size();

	LR_vector dra1 = _distance(p2a_vec, p1a_vec);
	LR_vector dra2 = _distance(p2a_vec, p3a_vec);

	LR_vector dra1_normed = dra1 / dra1.module();
	LR_vector dra2_normed = dra2 / dra2.module();

	number dot_product = dra1_normed * dra2_normed;
	number angle = std::acos(dot_product);
	int ix_left = std::floor((angle - this->xmin) / this->dX);
	int ix_right = ix_left + 1;

	number my_potential = 0;

	if((ix_left < 0) || (ix_right > N_grid - 1)) {
		std::cout << "off grid!" << std::endl;
	}

	else {
		my_potential = interpolatePotential(angle, dX, xmin, potential_grid);
	}

	number total_factor = _p1a_ptr.size() + _p2a_ptr.size() + _p3a_ptr.size();

	return my_potential / total_factor;
}

/*
 * MD_CPUBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MD_CPUBackend.h"
#include "Thermostats/ThermostatFactory.h"
#include "Thermostats/NoThermostat.h"
#include "MCMoves/MoveFactory.h"

MD_CPUBackend::MD_CPUBackend() :
				MDBackend() {
	_thermostat = nullptr;
	_V_move = nullptr;
	_stress_tensor_avg_every = -1;
	_stress_tensor_counter = 0;
}

MD_CPUBackend::~MD_CPUBackend() {

}

void MD_CPUBackend::get_settings(input_file &inp) {
	MDBackend::get_settings(inp);

	if(_use_barostat) {
		std::string move_name("volume");
		if(_barostat_molecular) {
			move_name = "molecule_volume";
		}
		auto str_inp = Utils::sformat("type = %s\ndelta = %lf\nisotropic = %d", move_name.c_str(), _delta_L, (int) _barostat_isotropic);
		input_file *move_inp = Utils::get_input_file_from_string(str_inp);
		_V_move = MoveFactory::make_move(*move_inp, inp);
		delete move_inp;
	}

	getInputBool(&inp, "MD_use_builtin_langevin_thermostat", &_use_builtin_langevin_thermostat, 0);
	if(_use_builtin_langevin_thermostat) {
		OX_LOG(Logger::LOG_INFO, "Using the built-in Langevin thermostat");
		_thermostat = std::shared_ptr<BaseThermostat>(new NoThermostat());

		number gamma;
		getInputNumber(&inp, "MD_langevin_gamma", &gamma, 1);
		_langevin_c1 = std::exp(-gamma * (_dt / 2.));
		_langevin_c2 = std::sqrt((1. - SQR(_langevin_c1)) * _T);
	}
	else {
		_thermostat = ThermostatFactory::make_thermostat(inp, _box.get());
		_thermostat->get_settings(inp);
	}
}

void MD_CPUBackend::init() {
	MDBackend::init();
	_thermostat->init();
	if(_use_barostat) {
		_V_move->init();
	}

	_compute_forces();
}

void MD_CPUBackend::_first_step() {
	std::vector<int> particles_with_warning;
	LR_vector dr;
	for(auto p : _particles) {
		if(_use_builtin_langevin_thermostat) {
			LR_vector p_plus = _langevin_c1 * p->vel + _langevin_c2 * LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian());
			LR_vector dv = p->force * (_dt * (number) 0.5);
			p->vel = p_plus + dv;
			dr = (p_plus + dv) * _dt;
		}
		else {
			p->vel += p->force * (_dt * (number) 0.5);
			dr = p->vel * _dt;
		}
		if(dr.norm() > 0.01) {
			particles_with_warning.push_back(p->index);
		}
		p->pos += dr;
		// if Lees-Edwards boundaries are enabled, we have to check for crossings along the y axis
		if(_lees_edwards) {
			const LR_vector &L = _box->box_sides();
			int y_new = floor(p->pos.y / L.y);
			int y_old = floor((p->pos.y - dr.y) / L.y);
			// we crossed the boundary along y
			if(y_new != y_old) {
				number delta_x = _shear_rate * L.y * current_step() * _dt;
				delta_x -= floor(delta_x / L.x) * L.x;
				if(y_new > y_old) {
					p->pos.x -= delta_x;
					p->pos.y -= L.y;
					p->vel.x -= _shear_rate * L.y;
				}
				else {
					p->pos.x += delta_x;
					p->pos.y += L.y;
					p->vel.x += _shear_rate * L.y;
				}
			}
		}

		if(p->is_rigid_body()) {
			p->L += p->torque * (_dt * (number) 0.5);
			// update of the orientation
			number norm = p->L.module();
			LR_vector LVersor(p->L / norm);

			number sintheta = sin(_dt * norm);
			number costheta = cos(_dt * norm);
			number olcos = 1. - costheta;

			number xyo = LVersor[0] * LVersor[1] * olcos;
			number xzo = LVersor[0] * LVersor[2] * olcos;
			number yzo = LVersor[1] * LVersor[2] * olcos;
			number xsin = LVersor[0] * sintheta;
			number ysin = LVersor[1] * sintheta;
			number zsin = LVersor[2] * sintheta;

			LR_matrix R(LVersor[0] * LVersor[0] * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, LVersor[1] * LVersor[1] * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, LVersor[2] * LVersor[2] * olcos + costheta);

			p->orientation = p->orientation * R;
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
		}

		p->set_initial_forces(current_step(), _box.get());

		_lists->single_update(p);
	}

	if(particles_with_warning.size() > 0) {
		std::stringstream ss;
		for(auto idx : particles_with_warning) {
			ss << idx << " ";
		}
		OX_LOG(Logger::LOG_WARNING, "The following particles had a displacement greater than 0.1 in this step: %s", ss.str().c_str());
	}
}

void MD_CPUBackend::_compute_forces() {
	_interaction->begin_energy_and_force_computation();

	_U = (number) 0;
	for(auto p : _particles) {
		for(auto &pair : p->affected) {
			if(pair.first == p) {
				_U += _interaction->pair_interaction_bonded(pair.first, pair.second, true, true);
			}
		}

		for(auto q : _lists->get_neigh_list(p)) {
			_U += _interaction->pair_interaction_nonbonded(p, q, true, true);
		}
	}
}

void MD_CPUBackend::_second_step() {
	for(auto p : _particles) {
		p->vel += p->force * _dt * (number) 0.5f;
		if(_use_builtin_langevin_thermostat) {
			p->vel = _langevin_c1 * p->vel + _langevin_c2 * LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian());
		}

		if(p->is_rigid_body()) {
			p->L += p->torque * _dt * (number) 0.5f;
		}
	}
}

void MD_CPUBackend::_update_backend_info() {

}

void MD_CPUBackend::sim_step() {
	_mytimer->resume();

	_timer_first_step->resume();
	_first_step();
	_timer_first_step->pause();

	_timer_lists->resume();
	if(!_lists->is_updated()) {
		_lists->global_update();
		_N_updates++;
	}
	_timer_lists->pause();

	if(_is_barostat_active()) {
		_timer_barostat->resume();
		_V_move->apply(current_step());
		_barostat_acceptance = _V_move->get_acceptance();
		_timer_barostat->pause();
	}

	_timer_forces->resume();
	_compute_forces();
	_second_step();

	_timer_forces->pause();

	_timer_thermostat->resume();
	_thermostat->apply(_particles, current_step());
	_timer_thermostat->pause();

	_mytimer->pause();
}

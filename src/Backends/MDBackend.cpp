/*
 * MDBackend.cpp
 *
 *  Created on: 23/nov/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MDBackend.h"
#include "../Observables/ObservableOutput.h"

MDBackend::MDBackend() :
				SimBackend() {
	_timer_first_step = _timer_forces = _timer_lists = _timer_thermostat = _timer_barostat = NULL;

	_lees_edwards = false;
	_shear_rate = -0.f;
}

MDBackend::~MDBackend() {

}

void MDBackend::get_settings(input_file &inp) {
	SimBackend::get_settings(inp);

	// refresh of initial velocities
	getInputBool(&inp, "refresh_vel", &_refresh_velocities, 0);
	getInputBool(&inp, "reset_initial_com_momentum", &_reset_initial_com_momentum, 0);
	getInputBool(&inp, "reset_com_momentum", &_reset_com_momentum, 0);
	if(_reset_com_momentum && !this->_enable_fix_diffusion) throw oxDNAException("reset_com_momentum requires fix_diffusion = true");

	getInputBool(&inp, "lees_edwards", &_lees_edwards, 0);
	if(_lees_edwards) {
		getInputNumber(&inp, "lees_edwards_shear_rate", &_shear_rate, 1);
		if(_shear_rate <= 0.) throw oxDNAException("lees_edwards_shear_rate should be > 0");
		OX_LOG(Logger::LOG_INFO, "Using Lees-Edwards boundary conditions with shear rate %lf", _shear_rate);
	}

	getInputNumber(&inp, "dt", &_dt, 1);

	getInputBool(&inp, "use_barostat", &_use_barostat, 0);
	if(_use_barostat) {
		if(_lees_edwards) {
			throw oxDNAException("Lees-Edwards boundaries are not compatible with isobaric simulations");
		}
		getInputBool(&inp, "barostat_isotropic", &_barostat_isotropic, 0);
		getInputNumber(&inp, "P", &this->_P, 1);
		getInputNumber(&inp, "delta_L", &_delta_L, 1);
		getInputNumber(&inp, "barostat_probability", &_barostat_probability, 1);
		getInputBool(&inp, "barostat_molecular", &_barostat_molecular, 0);
		std::string barostat_type("atomic");
		if(_barostat_molecular) {
			barostat_type = "molecular";
		}

		OX_LOG(Logger::LOG_INFO, "Enabling the MC-like %s barostat (isotropic = %d) with P = %lf, delta_L = %lf and activation probability = %lf", barostat_type.c_str(), _barostat_isotropic, this->_P, _delta_L, _barostat_probability);
	}

	char energy_file[512];
	llint print_every;
	getInputString(&inp, "energy_file", energy_file, 1);
	getInputLLInt(&inp, "print_energy_every", &print_every, 1);

	// we build the default stream of observables;
	// we use a helper string for that
	string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", energy_file, print_every);
	_obs_output_file = std::make_shared<ObservableOutput>(fake);
	add_output(_obs_output_file);
	_obs_output_file->add_observable("type = step\nunits = MD");
	_obs_output_file->add_observable("type = total_energy");
	if(_use_barostat) {
		_obs_output_file->add_observable("type = density");
	}
	_obs_output_file->add_observable("type = backend_info");

	// now we do the same thing for stdout
	int no_stdout_energy = 0;
	getInputBoolAsInt(&inp, "no_stdout_energy", &no_stdout_energy, 0);
	if(!no_stdout_energy) {
		fake = Utils::sformat("{\n\tname = stdout\n\tprint_every = %lld\n}\n", print_every);
		_obs_output_stdout = std::make_shared<ObservableOutput>(fake);
		add_output(_obs_output_stdout);
		_obs_output_stdout->add_observable("type = step");
		_obs_output_stdout->add_observable("type = step\nunits = MD");
		_obs_output_stdout->add_observable("type = total_energy");
		if(_use_barostat) {
			_obs_output_stdout->add_observable("type = density");
		}
		_obs_output_stdout->add_observable("type = backend_info");
	}
}

void MDBackend::init() {
	SimBackend::init();

	if(_refresh_velocities) _generate_vel();
	else {
		for(auto p: _particles) {
			if(p->L.module() < 1.e-10) {
				throw oxDNAException("Particle %i has 0 angular momentum in initial configuration.\n\tset \"refresh_vel = true\" in input file. Aborting now.", p->index);
			}
		}
	}

	if(_reset_initial_com_momentum) {
		OX_LOG(Logger::LOG_INFO, "Setting the centre of mass' momentum to 0");
		_reset_momentum();
	}

	_timer_first_step = TimingManager::instance()->new_timer(string("First Step"), string("SimBackend"));
	_timer_forces = TimingManager::instance()->new_timer(string("Forces"), string("SimBackend"));
	_timer_thermostat = TimingManager::instance()->new_timer(string("Thermostat"), string("SimBackend"));
	_timer_lists = TimingManager::instance()->new_timer(string("Lists"), string("SimBackend"));
	if(_use_barostat) {
		_timer_barostat = TimingManager::instance()->new_timer(string("Barostat"), string("SimBackend"));
	}
}

bool MDBackend::_is_barostat_active() {
	if(!_use_barostat) {
		return false;
	}
	return _barostat_probability > drand48();
}

void MDBackend::_reset_momentum() {
	LR_vector com_v(0, 0, 0);
	for(auto p: _particles) {
		com_v += p->vel;
	}
	com_v /= N();

	for(auto p: _particles) {
		p->vel -= com_v;
	}
}

void MDBackend::_generate_vel() {
	OX_LOG(Logger::LOG_INFO, "Using randomly distributed velocities");

	number rescale_factor = sqrt(this->_T);
	number initial_K = 0;
	for(auto p: _particles) {

		p->vel.x = Utils::gaussian() * rescale_factor;
		p->vel.y = Utils::gaussian() * rescale_factor;
		p->vel.z = Utils::gaussian() * rescale_factor;

		p->L.x = Utils::gaussian() * rescale_factor;
		p->L.y = Utils::gaussian() * rescale_factor;
		p->L.z = Utils::gaussian() * rescale_factor;

		if(_lees_edwards) {
			number Ly = this->_box->box_sides().y;
			number y_in_box = p->pos.y - floor(p->pos.y/Ly)*Ly - 0.5*Ly;
			p->vel.x += y_in_box*_shear_rate;
		}

		initial_K += (p->vel.norm() + p->L.norm()) * 0.5;
	}

	OX_LOG(Logger::LOG_INFO, "Initial kinetic energy: %f", initial_K);
}

void MDBackend::fix_diffusion() {
	if(_reset_com_momentum) MDBackend::_reset_momentum();
	SimBackend::fix_diffusion();
}

void MDBackend::print_observables() {
	if(_use_barostat) {
		this->_backend_info.insert(0, Utils::sformat(" %5.3lf", _barostat_acceptance));
	}

	SimBackend::print_observables();
}

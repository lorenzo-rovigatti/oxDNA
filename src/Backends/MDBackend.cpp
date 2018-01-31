/*
 * MDBackend.cpp
 *
 *  Created on: 23/nov/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MDBackend.h"
#include "../Observables/ObservableOutput.h"

template<typename number>
MDBackend<number>::MDBackend() : SimBackend<number>(), _refresh_velocities(false) {
	this->_sim_type = SIM_MD;
	_reset_initial_com_momentum = false;
	_reset_com_momentum = false;
	_use_barostat = false;
	_barostat_isotropic = true;
	_barostat_acceptance = 0.;

	_timer_first_step = _timer_forces = _timer_lists = _timer_thermostat = _timer_barostat = NULL;

	_lees_edwards = false;
	_shear_rate = -0.f;
}

template<typename number>
MDBackend<number>::~MDBackend() {

}

template<typename number>
void MDBackend<number>::get_settings(input_file &inp) {
	SimBackend<number>::get_settings(inp);

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

	getInputNumber<number>(&inp, "dt", &_dt, 1);

	getInputBool(&inp, "use_barostat", &_use_barostat, 0);
	if(_use_barostat) {
		if(_lees_edwards) throw oxDNAException("Lees-Edwards boundaries are not compatible with isobaric simulations");
		getInputBool(&inp, "barostat_isotropic", &_barostat_isotropic, 0);
		getInputNumber<number>(&inp, "P", &this->_P, 1);
		getInputNumber<number>(&inp, "delta_L", &_delta_L, 1);
		_barostat_probability = 1.;
		getInputNumber<number>(&inp, "barostat_probability", &_barostat_probability, 0);
		OX_LOG(Logger::LOG_INFO, "Enabling the MC-like barostat (isotropic = %d) with P = %lf, delta_L = %lf and activation probability = %lf", _barostat_isotropic, this->_P, _delta_L, _barostat_probability);
	}

	char energy_file[512];
	llint print_every;
	getInputString(&inp, "energy_file", energy_file, 1);
	getInputLLInt(&inp, "print_energy_every", &print_every, 1);

	// we build the default stream of observables;
	// we use a helper string for that
	string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", energy_file, print_every);
	this->_obs_output_file = new ObservableOutput<number>(fake, inp);
	this->_obs_outputs.push_back(this->_obs_output_file);
	this->_obs_output_file->add_observable("type = step\nunits = MD");
	this->_obs_output_file->add_observable("type = total_energy");
	if(_use_barostat) this->_obs_output_file->add_observable("type = density");
	this->_obs_output_file->add_observable("type = backend_info");

	// now we do the same thing for stdout
	int no_stdout_energy = 0;
	getInputBoolAsInt(&inp, "no_stdout_energy", &no_stdout_energy, 0);
	if(!no_stdout_energy) {
		fake = Utils::sformat("{\n\tname = stdout\n\tprint_every = %lld\n}\n", print_every);
		this->_obs_output_stdout = new ObservableOutput<number>(fake, inp);
		this->_obs_outputs.push_back(this->_obs_output_stdout);
		this->_obs_output_stdout->add_observable("type = step");
		this->_obs_output_stdout->add_observable("type = step\nunits = MD");
		this->_obs_output_stdout->add_observable("type = total_energy");
		if(_use_barostat) this->_obs_output_stdout->add_observable("type = density");
		this->_obs_output_stdout->add_observable("type = backend_info");
	}
}

template<typename number>
void MDBackend<number>::init() {
	SimBackend<number>::init();

	if(_refresh_velocities) _generate_vel();
	else {
	    for(int i = 0; i < this->_N; i ++) {
	    	if(this->_particles[i]->L.module() < 1.e-10)
	    		throw oxDNAException("Particle %i has 0 angular momentum in initial configuration.\n\tset \"refresh_vel = 1\" in input file. Aborting now.", i);
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
	if(_use_barostat) _timer_barostat = TimingManager::instance()->new_timer(string("Barostat"), string("SimBackend"));
}

template<typename number>
bool MDBackend<number>::_is_barostat_active() {
	if(!_use_barostat) return false;
	return _barostat_probability > drand48();
}

template<typename number>
void MDBackend<number>::_reset_momentum() {
	LR_vector<number> com_v(0, 0, 0);
	for(int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		com_v += p->vel;
	}
	com_v /= this->_N;

	for(int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		p->vel -= com_v;
	}
}

template<typename number>
void MDBackend<number>::_generate_vel() {
	OX_LOG(Logger::LOG_INFO, "Using randomly distributed velocities");

	number rescale_factor = sqrt(this->_T);
	number initial_K = 0;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->vel.x = Utils::gaussian<number>() * rescale_factor;
		p->vel.y = Utils::gaussian<number>() * rescale_factor;
		p->vel.z = Utils::gaussian<number>() * rescale_factor;

		p->L.x = Utils::gaussian<number>() * rescale_factor;
		p->L.y = Utils::gaussian<number>() * rescale_factor;
		p->L.z = Utils::gaussian<number>() * rescale_factor;

		if(_lees_edwards) {
			number Ly = this->_box->box_sides().y;
			number y_in_box = p->pos.y - floor(p->pos.y/Ly)*Ly - 0.5*Ly;
			p->vel.x += y_in_box*_shear_rate;
		}

		initial_K += (p->vel.norm() + p->L.norm()) * 0.5;
	}

	OX_LOG(Logger::LOG_INFO, "Initial kinetic energy: %f", initial_K);
}

template<typename number>
void MDBackend<number>::fix_diffusion() {
	if(_reset_com_momentum) MDBackend<number>::_reset_momentum();
	SimBackend<number>::fix_diffusion();
}

template<typename number>
void MDBackend<number>::print_observables(llint curr_step) {
	if(_use_barostat) this->_backend_info.insert(0, Utils::sformat(" %5.3lf", _barostat_acceptance));

	SimBackend<number>::print_observables(curr_step);
}

template class MDBackend<float>;
template class MDBackend<double>;

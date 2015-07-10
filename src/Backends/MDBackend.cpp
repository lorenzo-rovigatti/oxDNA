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

	getInputNumber<number>(&inp, "dt", &_dt, 1);

	char energy_file[512];
	llint print_every;
	getInputString(&inp, "energy_file", energy_file, 1);
	getInputLLInt(&inp, "print_energy_every", &print_every, 1);

	// we build the default stream of observables;
	// we build a helper string for that
	std::string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", energy_file, print_every);
	this->_obs_output_file = new ObservableOutput<number>(fake, inp);
	this->_obs_outputs.push_back(this->_obs_output_file);
	this->_obs_output_file->add_observable("type = step\nunits = MD");
	this->_obs_output_file->add_observable("type = total_energy");
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

	_timer_first_step = TimingManager::instance()->new_timer(std::string("First Step"), std::string("SimBackend"));
	_timer_forces = TimingManager::instance()->new_timer(std::string("Forces"), std::string("SimBackend"));
	_timer_thermostat = TimingManager::instance()->new_timer(std::string("Thermostat"), std::string("SimBackend"));
	_timer_lists = TimingManager::instance()->new_timer(std::string("Lists"), std::string("SimBackend"));
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

		initial_K += (p->vel.norm() + p->L.norm()) * 0.5;
	}

	OX_LOG(Logger::LOG_INFO, "Initial kinetic energy: %f", initial_K);
}

template<typename number>
void MDBackend<number>::fix_diffusion() {
	if(_reset_com_momentum) MDBackend<number>::_reset_momentum();
	SimBackend<number>::fix_diffusion();
}

template class MDBackend<float>;
template class MDBackend<double>;

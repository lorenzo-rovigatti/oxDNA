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

	// initialize the messages for the timings output
	this->_timer_msgs_number = 6;
	strncpy(this->_timer_msgs[0], "MD step", 256);
	strncpy(this->_timer_msgs[1], "First step", 256);
	strncpy(this->_timer_msgs[2], "Hilbert sorting", 256);
	strncpy(this->_timer_msgs[3], "Lists update", 256);
	strncpy(this->_timer_msgs[4], "Forces + second step", 256);
	strncpy(this->_timer_msgs[5], "Thermostat", 256);
}

template<typename number>
MDBackend<number>::~MDBackend() {

}

template<typename number>
void MDBackend<number>::get_settings(input_file &inp) {
	SimBackend<number>::get_settings(inp);

	// refresh of initial velocities
	int resc;
	if(getInputBoolAsInt(&inp, "refresh_vel", &resc, 0) == KEY_FOUND) _refresh_velocities = resc;

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
void MDBackend<number>::init(char conf_filename[256]) {
	SimBackend<number>::init(conf_filename);

	if(_refresh_velocities) _generate_vel();
	else {
	    for(int i = 0; i < this->_N; i ++) {
	    	if(this->_particles[i]->L.module() < 1.e-10)
	    		throw oxDNAException("Particle %i has 0 angular momentum in initial configuration.\n\tset \"refresh_vel = 1\" in input file. Aborting now.", i);
		}
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

template class MDBackend<float>;
template class MDBackend<double>;

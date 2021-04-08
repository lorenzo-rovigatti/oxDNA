/*
 * MDBackend.cpp
 *
 *  Created on: 23/nov/2010
 *      Author: lorenzo
 */

#include "MDBackend.h"
#include "IOManager.h"

template<typename number>
MDBackend<number>::MDBackend(IOManager *IO) : SimBackend<number>(IO), _refresh_velocities(false), _newtonian_steps(0), _thermostat(THERMOSTAT_NO), _pt(0), _diff_coeff(0){
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

template<>
void MDBackend<float>::_get_number_settings(input_file &inp) {
	getInputFloat(&inp, "dt", &_dt, 1);
	if(_thermostat == this->THERMOSTAT_JOHN) {
		if(getInputFloat(&inp, "pt", &_pt, 0) == KEY_NOT_FOUND) {
			if(getInputFloat(&inp, "diff_coeff", &_diff_coeff, 0) == KEY_NOT_FOUND)
				this->_IO->die("pt or diff_coeff must be specified for the John thermostat");
		}
	}
}

template<>
void MDBackend<double>::_get_number_settings(input_file &inp) {
	getInputDouble(&inp, "dt", &_dt, 1);
	if(_thermostat == this->THERMOSTAT_JOHN) {
		if(getInputDouble(&inp, "pt", &_pt, 0) == KEY_NOT_FOUND) {
			if(getInputDouble(&inp, "diff_coeff", &_diff_coeff, 0) == KEY_NOT_FOUND)
				this->_IO->die("pt or diff_coeff must be specified for the John thermostat");
		}
	}
}

template<typename number>
void MDBackend<number>::get_settings(input_file &inp) {
	char key[256];

	SimBackend<number>::get_settings(inp);

	int resc;
	if(getInputInt(&inp, "refresh_vel", &resc, 0) == KEY_FOUND) _refresh_velocities = (bool) resc;

	if(getInputString(&inp, "thermostat", key, 0) == KEY_FOUND) {
		if(!strcmp(key, "no")) _thermostat = this->THERMOSTAT_NO;
		else {
			if(!strcmp(key, "john")) _thermostat = this->THERMOSTAT_JOHN;
			else if(!strcmp(key, "refresh")) _thermostat = this->THERMOSTAT_REFRESH;
			else this->_IO->die("Thermostat '%s' not supported", key);

			getInputInt(&inp, "newtonian_steps", &_newtonian_steps, 1);
			if(_newtonian_steps < 1) this->_IO->die("'newtonian_steps' must be > 0");
		}
	}

	_get_number_settings(inp);
}

template<typename number>
//void MDBackend<number>::init(ifstream &conf_input) {
void MDBackend<number>::init(char conf_filename[256]) {
	//SimBackend<number>::init(conf_input);
	SimBackend<number>::init(conf_filename);

	if(_thermostat == this->THERMOSTAT_JOHN) {
		if(_pt == (number) 0.) _pt = (2 * this->_T *  _newtonian_steps * _dt)/(this->_T * _newtonian_steps * _dt + 2 * _diff_coeff);
		if(_pt >= (number) 1.) this->_IO->die("pt (%f) must be smaller than 1", _pt);

		// initialize pr (considering Dr = 3Dt)
		_diff_coeff = this->_T * _newtonian_steps * _dt * (1./_pt - 1./2.);
		_pr = (2 * this->_T *  _newtonian_steps * _dt)/(this->_T * _newtonian_steps * _dt + 2 * 3 * _diff_coeff);
	}

	_rescale_factor = sqrt(this->_T);

	if(_refresh_velocities) _generate_vel();
	else {
	    for (int i = 0; i < this->_N; i ++) {
		if (this->_particles[i].L.module() < 1.e-10) this->_IO->die("Particle %i has 0 angular momentum in initial configuration.\n\tset \"refresh_vel = 1\" in input file. Aborting now.", i);
		}
	}

	if(_thermostat == this->THERMOSTAT_JOHN) this->_IO->log(this->_IO->LOG_INFO, "John's pt: %f, pr: %lf", _pt, _pr);
}

template<typename number>
void MDBackend<number>::_generate_vel() {
	this->_IO->log(this->_IO->LOG_INFO, "Using randomly distributed velocities");

	number initial_K = 0;
	for(int i = 0; i < this->_N; i++) {
		this->_particles[i].vel.x = Utils::gaussian<number>() * this->_rescale_factor;
		this->_particles[i].vel.y = Utils::gaussian<number>() * this->_rescale_factor;
		this->_particles[i].vel.z = Utils::gaussian<number>() * this->_rescale_factor;

		this->_particles[i].L.x = Utils::gaussian<number>() * this->_rescale_factor;
		this->_particles[i].L.y = Utils::gaussian<number>() * this->_rescale_factor;
		this->_particles[i].L.z = Utils::gaussian<number>() * this->_rescale_factor;

		initial_K += (this->_particles[i].vel.norm() + this->_particles[i].L.norm()) * 0.5;
	}

	this->_IO->log(this->_IO->LOG_INFO, "Initial kinetic energy: %f", initial_K);
}

template<typename number>
void MDBackend<number>::print_energy(llint curr_step) {
	this->_IO->print_energy(*this, curr_step);
}

template<typename number>
void MDBackend<number>::print_conf(llint curr_step, bool reduced, bool only_last) {
	if(reduced == false) this->_IO->print_conf(*this, curr_step, only_last);
	else this->_IO->print_reduced_conf(*this, curr_step);
}

template class MDBackend<float>;
template class MDBackend<double>;

/*
 * MD_CPUBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MD_CPUBackend.h"
#include "./Thermostats/ThermostatFactory.h"

template<typename number>
MD_CPUBackend<number>::MD_CPUBackend() : MDBackend<number>() {
	this->_is_CUDA_sim = false;
	_thermostat = NULL;
}

template<typename number>
MD_CPUBackend<number>::~MD_CPUBackend() {
	delete _thermostat;
}

template<typename number>
void MD_CPUBackend<number>::_first_step(llint cur_step) {
	bool is_warning = false;
	std::vector<int> w_ps;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->vel += p->force * this->_dt * (number) 0.5;
		LR_vector<number> dr = p->vel * this->_dt;
		if(dr.norm() > 0.01) {
			is_warning = true;
			w_ps.push_back(p->index);
		}
		p->pos += dr;

		if(p->is_rigid_body()) {
			p->L += p->torque * this->_dt * (number) 0.5;
			// update of the orientation
			number norm = p->L.module();
			LR_vector<number> LVersor(p->L / norm);

			number sintheta = sin(this->_dt * norm);
			number costheta = cos(this->_dt * norm);
			number olcos = 1. - costheta;

			number xyo = LVersor[0] * LVersor[1] * olcos;
			number xzo = LVersor[0] * LVersor[2] * olcos;
			number yzo = LVersor[1] * LVersor[2] * olcos;
			number xsin = LVersor[0] * sintheta;
			number ysin = LVersor[1] * sintheta;
			number zsin = LVersor[2] * sintheta;

			LR_matrix<number> R(LVersor[0] * LVersor[0] * olcos + costheta, xyo - zsin, xzo + ysin,
						xyo + zsin, LVersor[1] * LVersor[1] * olcos + costheta, yzo - xsin,
						xzo - ysin, yzo + xsin, LVersor[2] * LVersor[2] * olcos + costheta);

			p->orientation = p->orientation * R;
			// set back, base and stack positions
			p->set_positions();
			p->orientationT = p->orientation.get_transpose();
			p->torque = LR_vector<number>((number) 0, (number) 0, (number) 0);
		}

		p->set_initial_forces(cur_step, this->_box_side);

		this->_lists->single_update(p);
	}

	if(is_warning) {
		std::stringstream ss;
		for(vector<int>::iterator it = w_ps.begin(); it != w_ps.end(); it++) ss << *it << " ";
		OX_LOG(Logger::LOG_WARNING, "The following particles had a displacement greater than 0.1 in this step: %s", ss.str().c_str());
	}
}

template<typename number>
void MD_CPUBackend<number>::_compute_forces() {
	this->_U = this->_U_hydr = (number) 0;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		this->_U += this->_interaction->pair_interaction_bonded(p, P_VIRTUAL, NULL, true);

		std::vector<BaseParticle<number> *> neighs = this->_lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle<number> *q = neighs[n];
			this->_U += this->_interaction->pair_interaction_nonbonded(p, q, NULL, true);
		}
	}
}

template<typename number>
void MD_CPUBackend<number>::_second_step() {
	this->_K = (number) 0.f;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->vel += p->force * this->_dt * (number) 0.5f;
		if(p->is_rigid_body()) p->L += p->torque * this->_dt * (number) 0.5f;

		this->_K += (p->vel.norm() + p->L.norm()) * (number) 0.5f;
	}
}

template<typename number>
void MD_CPUBackend<number>::sim_step(llint curr_step) {
	this->_mytimer->resume();

	this->_timer_first_step->resume();
	_first_step(curr_step);
	this->_timer_first_step->pause();

	this->_timer_lists->resume();
	if(!this->_lists->is_updated()) {
		this->_lists->global_update();
		this->_N_updates++;
	}
	this->_timer_lists->pause();

	this->_timer_forces->resume();
	_compute_forces();
	_second_step();
	this->_timer_forces->pause();

	this->_timer_thermostat->resume();
	_thermostat->apply(this->_particles, curr_step);
	this->_timer_thermostat->pause();

	this->_mytimer->pause();
}
template<typename number>
void MD_CPUBackend<number>::get_settings (input_file &inp) {
	MDBackend<number>::get_settings(inp);
	_thermostat = ThermostatFactory::make_thermostat<number>(inp, this->_box_side);
	_thermostat->get_settings(inp);
}

template<typename number>
void MD_CPUBackend<number>::init() {
	MDBackend<number>::init();
	_thermostat->init (this->_N);

	_compute_forces();
}

template class MD_CPUBackend<float>;
template class MD_CPUBackend<double>;

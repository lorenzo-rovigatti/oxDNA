/*
 * MinBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MinBackend.h"
#include "./Thermostats/ThermostatFactory.h"

template<typename number>
MinBackend<number>::MinBackend() : MDBackend<number>() {
	this->_is_CUDA_sim = false;

	_max_step = 0.005;
}

template<typename number>
MinBackend<number>::~MinBackend() {

}

template<typename number>
void MinBackend<number>::get_settings (input_file &inp) {
	MDBackend<number>::get_settings(inp);

	getInputNumber (&inp, "minimization_max_step", &_max_step, 0);

}

template<typename number>
void MinBackend<number>::init () {
	MDBackend<number>::init();

	_compute_forces();
	
	OX_LOG(Logger::LOG_INFO, "Using max step for minimization = %g", _max_step);
}


template<typename number>
void MinBackend<number>::_compute_forces() {
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
void MinBackend<number>::_evolve () {
	// find maximum force
	number max_f = -1.f;
	number max_t = -1.f;

	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		number tmp = p->force.norm();
		if (tmp > max_f) max_f = tmp;
		
		tmp = p->torque.norm();
		if (tmp > max_t) max_t = tmp;

	}

	// now let's look for the appropriate factor;
	number fact_r = 0.f;
	number fact_l = 0.f;
	fact_r = sqrt(max_f) / _max_step;
	fact_l = sqrt(max_t) / _max_step;

	//if (fact_r < 1.f) fact_r = 1.f;
	//if (fact_l < 1.f) fact_l = 1.f;

	// we evolve all the particles' position
	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		p->pos = p->pos + p->force / fact_r;
		if (p->is_rigid_body()) {
			// update of the orientation
			number norm = p->torque.module();
			LR_vector<number> LVersor(p->torque / norm);

			number sintheta = sin(norm / fact_l);
			number costheta = cos(norm / fact_l);
			number olcos = 1.f - costheta;

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
		this->_lists->single_update(p);
	}

	return;
}

template<typename number>
void MinBackend<number>::sim_step(llint curr_step) {
	this->_mytimer->resume();
	
	for(int i = 0; i < this->_N; i++) this->_particles[i]->set_initial_forces(curr_step, this->_box_side);

	this->_timer_lists->resume();
	if(!this->_lists->is_updated()) {
		this->_lists->global_update();
		this->_N_updates++;
	}
	this->_timer_lists->pause();

	this->_timer_forces->resume();
	_compute_forces();
	this->_timer_forces->pause();
	
	// here we normalize forces
	//for (int i = 0; i < this->_N; i ++) {
	//	BaseParticle<number> * p = this->_particles[i];
	//	if (p->force.module() > 0.1) p->force = p->force * (0.1 / p->force.module());
	//}
	
	_evolve();

	this->_mytimer->pause();
}

template class MinBackend<float>;
template class MinBackend<double>;


/*
 * MinBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include "MinBackend.h"

#include "../Observables/ObservableOutput.h"
#include "./Thermostats/ThermostatFactory.h"

#include <sstream>

MinBackend::MinBackend() :
				MDBackend() {
	_max_step = 0.005;
}

MinBackend::~MinBackend() {

}

void MinBackend::get_settings(input_file &inp) {
	MDBackend::get_settings(inp);

	getInputNumber(&inp, "minimization_max_step", &_max_step, 0);

}

void MinBackend::init() {
	MDBackend::init();

	_compute_forces();

	OX_LOG(Logger::LOG_INFO, "Using max step for minimization = %g", _max_step);
}

void MinBackend::_compute_forces() {
	_interaction->begin_energy_and_force_computation();

	_U = (number) 0;
	for(auto p: _particles) {
		for(auto &pair : p->affected) {
			if(pair.first == p) {
				_U += _interaction->pair_interaction_bonded(pair.first, pair.second, true, true);
			}
		}

		std::vector<BaseParticle *> neighs = _lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			_U += _interaction->pair_interaction_nonbonded(p, q, true, true);
		}
	}
}

void MinBackend::_evolve() {
	// find maximum force
	number max_f = -1.f;
	number max_t = -1.f;

	for(auto p: _particles) {
		number tmp = p->force.norm();
		if(tmp > max_f) max_f = tmp;

		tmp = p->torque.norm();
		if(tmp > max_t) max_t = tmp;

	}

	// now let's look for the appropriate factor;
	number fact_r = 0.f;
	number fact_l = 0.f;
	fact_r = sqrt(max_f) / _max_step;
	fact_l = sqrt(max_t) / _max_step;

	//if (fact_r < 1.f) fact_r = 1.f;
	//if (fact_l < 1.f) fact_l = 1.f;

	// we evolve all the particles' position
	for(auto p: _particles) {
		p->pos = p->pos + p->force / fact_r;
		if(p->is_rigid_body()) {
			// update of the orientation
			number norm = p->torque.module();
			LR_vector LVersor(p->torque / norm);

			number sintheta = sin(norm / fact_l);
			number costheta = cos(norm / fact_l);
			number olcos = 1.f - costheta;

			number xyo = LVersor[0] * LVersor[1] * olcos;
			number xzo = LVersor[0] * LVersor[2] * olcos;
			number yzo = LVersor[1] * LVersor[2] * olcos;
			number xsin = LVersor[0] * sintheta;
			number ysin = LVersor[1] * sintheta;
			number zsin = LVersor[2] * sintheta;

			LR_matrix R(LVersor[0] * LVersor[0] * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, LVersor[1] * LVersor[1] * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, LVersor[2] * LVersor[2] * olcos + costheta);

			p->orientation = p->orientation * R;
			// set back, base and stack positions
			p->set_positions();
			p->orientationT = p->orientation.get_transpose();
			p->torque = LR_vector((number) 0, (number) 0, (number) 0);
		}
		_lists->single_update(p);
	}

	return;
}

void MinBackend::sim_step() {
	_mytimer->resume();

	for(auto p: _particles) {
		p->set_initial_forces(current_step(), _box.get());
	}

	_timer_lists->resume();
	if(!_lists->is_updated()) {
		_lists->global_update();
		_N_updates++;
	}
	_timer_lists->pause();

	_timer_forces->resume();
	_compute_forces();
	_timer_forces->pause();

	// here we normalize forces
	//for (int i = 0; i < _N; i ++) {
	//	BaseParticle * p = _particles[i];
	//	if (p->force.module() > 0.1) p->force = p->force * (0.1 / p->force.module());
	//}

	_evolve();

	_mytimer->pause();
}


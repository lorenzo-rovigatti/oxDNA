/*
 * FIREBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */
#include "FIREBackend.h"

#include "../Observables/ObservableOutput.h"
#include "./Thermostats/ThermostatFactory.h"

#include <sstream>

FIREBackend::FIREBackend() :
				MDBackend() {

	_N_min = 5;
	_f_inc = 1.1;
	_f_dec = 0.5;
	_alpha_start = 0.1;
	_f_alpha = 0.99;
	_max_step = 0.005;

	_step_last_reset = 0;

	_dt_max = 0.005;
	_dt_restart = 1.e-16;

	_allow_dt_increase = false;
}

FIREBackend::~FIREBackend() {

}

void FIREBackend::get_settings(input_file &inp) {
	MDBackend::get_settings(inp);

	getInputNumber(&inp, "dt", &(_dt), 1);
	_dt_restart = _dt;

	getInputNumber(&inp, "minimization_max_step", &_max_step, 0);

}

void FIREBackend::init() {
	MDBackend::init();

	_compute_forces();

	_alpha = _alpha_start;
	_step_last_reset = 0;

	OX_LOG(Logger::LOG_INFO, "FIRE: Using max step for minimization = %g", _max_step);
}

void FIREBackend::_first_step() {
	// Spalletti
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

	if(fact_r < 1. && fact_l < 1.) _allow_dt_increase = true;
	else _allow_dt_increase = true;

	if(fact_r < 1.f) fact_r = 1.f;
	if(fact_l < 1.f) fact_l = 1.f;

	bool is_warning = false;
	std::vector<int> w_ps;
	std::vector<LR_vector> drs;
	std::vector<LR_matrix> dRs;
	std::vector<number> dthetas;

	for(auto p: _particles) {
		p->vel += p->force * _dt * (number) 0.5;
		LR_vector dr = p->vel * _dt;
		if(dr.norm() > 0.01) {
			is_warning = true;
			w_ps.push_back(p->index);
		}
		//p->pos += dr;
		drs.push_back(LR_vector(dr.x, dr.y, dr.z));

		if(p->is_rigid_body()) {
			p->L += p->torque * _dt * (number) 0.5;
			// update of the orientation
			number norm = p->L.module();
			LR_vector LVersor(p->L / norm);
			dthetas.push_back(_dt * norm);

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
			dRs.push_back(LR_matrix(LVersor[0] * LVersor[0] * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, LVersor[1] * LVersor[1] * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, LVersor[2] * LVersor[2] * olcos + costheta));

			//p->orientation = p->orientation * R;
			//p->orientationT = p->orientation.get_transpose();
			//p->set_positions();
			//p->torque = LR_vector((number) 0, (number) 0, (number) 0);
		}
	}

	_allow_dt_increase = true;
	number fact = -1.f;
	for(auto p: _particles) {
		LR_vector mydr = drs.back();
		drs.pop_back();
		if(mydr.module() > _max_step) {
			_allow_dt_increase = false;
			if((mydr.module() / _max_step) > fact) fact = mydr.module() / _max_step;
		}

		if(p->is_rigid_body()) {
			number mydtheta = dthetas.back();
			dthetas.pop_back();
			if(fabs(mydtheta) > _max_step) {
				_allow_dt_increase = false;
				if((fabs(mydtheta) / _max_step) > fact) fact = fabs(mydtheta) / _max_step;
			}
		}
	}

	if(fact > 1.f && _allow_dt_increase) throw oxDNAException("should never happen");

	if(fact < 1.f) fact = 1.f;
	printf("halfway fact: %g\n", fact);

	// here we increment coordinates
	for(auto p: _particles) {
		p->pos += p->vel * (_dt / fact);

		if(p->is_rigid_body()) {
			LR_matrix myR = dRs.back();
			dRs.pop_back();
			p->orientation = (p->orientation * myR) / fact;
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
			p->torque = LR_vector((number) 0, (number) 0, (number) 0);

		}
		p->set_initial_forces(current_step(), _box.get());

		_lists->single_update(p);
	}

	if(is_warning) {
		std::stringstream ss;
		for(vector<int>::iterator it = w_ps.begin(); it != w_ps.end(); it++)
			ss << *it << " ";
		OX_LOG(Logger::LOG_WARNING, "The following particles had a displacement greater than 0.1 in this step: %s", ss.str().c_str());
	}
}

void FIREBackend::_second_step() {
	for(auto p: _particles) {
		p->vel += p->force * _dt * (number) 0.5f;
		if(p->is_rigid_body()) {
			p->L += p->torque * _dt * (number) 0.5f;
		}
	}
}

void FIREBackend::_compute_forces() {
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

void FIREBackend::_evolve() {
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

	if(fact_r < 1.f) fact_r = 1.f;
	if(fact_l < 1.f) fact_l = 1.f;

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

void FIREBackend::sim_step() {
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

	_timer_forces->resume();
	_compute_forces();
	_second_step();
	_timer_forces->pause();

	_timer_thermostat->resume();

	// calculate P
	number P = 0.;
	for(auto p: _particles) {
		P += p->force * p->vel;
		if(p->is_rigid_body()) P += p->L * p->torque;
	}

	for(auto p: _particles) {
		p->vel = p->vel * (1. - _alpha) + p->force * (_alpha * p->vel.module() / p->force.module());
		if(p->is_rigid_body()) {
			p->L = p->L * (1. - _alpha) + p->torque * (_alpha * p->L.module() / p->torque.module());
		}
	}

	_step_last_reset += 1;

	if(P > (number) 0.f && _step_last_reset > _N_min && _allow_dt_increase) {
		_dt = _dt * _f_inc;
		if(_dt > _dt_max) _dt = _dt_max;
		_alpha = _alpha * _f_alpha;
	}
	else if(P <= (number) 0.f) {
		_dt = _dt * _f_dec;
		_alpha = _alpha_start;
		for(auto p: _particles) {
			p->vel.x = (number) 0.f;
			p->vel.y = (number) 0.f;
			p->vel.z = (number) 0.f;
			if(p->is_rigid_body()) {
				p->L.x = (number) 0.f;
				p->L.y = (number) 0.f;
				p->L.z = (number) 0.f;
			}
		}
		_step_last_reset = 0;
	}
	_timer_thermostat->pause();

	_mytimer->pause();
	printf("P, alpha, dt = %g, %g, %g\n", P, _alpha, _dt);
}

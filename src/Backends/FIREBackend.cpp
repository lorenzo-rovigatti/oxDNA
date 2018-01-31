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

template<typename number>
FIREBackend<number>::FIREBackend() : MDBackend<number>() {
	this->_is_CUDA_sim = false;

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

template<typename number>
FIREBackend<number>::~FIREBackend() {

}

template<typename number>
void FIREBackend<number>::get_settings (input_file &inp) {
	MDBackend<number>::get_settings(inp);

	getInputNumber (&inp, "dt", &(this->_dt), 1);
	_dt_restart = this->_dt;

	getInputNumber (&inp, "minimization_max_step", &_max_step, 0);

}

template<typename number>
void FIREBackend<number>::init () {
	MDBackend<number>::init();

	_compute_forces();

	_alpha = _alpha_start;
	_step_last_reset = 0;
	
	OX_LOG(Logger::LOG_INFO, "FIRE: Using max step for minimization = %g", _max_step);
}

template<typename number>
void FIREBackend<number>::_first_step(llint cur_step) {
	// Spalletti
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

	if (fact_r < 1. && fact_l < 1.) _allow_dt_increase = true;
	else _allow_dt_increase = true;

	if (fact_r < 1.f) fact_r = 1.f;
	if (fact_l < 1.f) fact_l = 1.f;

	bool is_warning = false;
	std::vector<int> w_ps;
	std::vector<LR_vector<number> > drs; 
	std::vector<LR_matrix<number> > dRs; 
	std::vector<number> dthetas;

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->vel += p->force * this->_dt * (number) 0.5;
		LR_vector<number> dr = p->vel * this->_dt;
		if(dr.norm() > 0.01) {
			is_warning = true;
			w_ps.push_back(p->index);
		}
		//p->pos += dr;
		drs.push_back(LR_vector<number> (dr.x, dr.y, dr.z));

		if(p->is_rigid_body()) {
			p->L += p->torque * this->_dt * (number) 0.5;
			// update of the orientation
			number norm = p->L.module();
			LR_vector<number> LVersor(p->L / norm);
			dthetas.push_back (this->_dt * norm);

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
			dRs.push_back (LR_matrix<number> (LVersor[0] * LVersor[0] * olcos + costheta, xyo - zsin, xzo + ysin,
						xyo + zsin, LVersor[1] * LVersor[1] * olcos + costheta, yzo - xsin,
						xzo - ysin, yzo + xsin, LVersor[2] * LVersor[2] * olcos + costheta));

			//p->orientation = p->orientation * R;
			//p->orientationT = p->orientation.get_transpose();
			//p->set_positions();
			//p->torque = LR_vector<number>((number) 0, (number) 0, (number) 0);
		}
	}

	_allow_dt_increase = true;
	number fact = -1.f;
	for (int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		LR_vector<number> mydr = drs.back();
		drs.pop_back();
		if (mydr.module() > _max_step) {
			_allow_dt_increase = false;
			if ((mydr.module() / _max_step) > fact)
				fact = mydr.module() / _max_step;
		}

		if (p->is_rigid_body()) {
			number mydtheta = dthetas.back();
			dthetas.pop_back();
			if (fabs(mydtheta) > _max_step) {
				_allow_dt_increase = false;
				if ((fabs(mydtheta) / _max_step) > fact)
					fact = fabs(mydtheta) / _max_step;
			}
		}
	}

	if (fact > 1.f && _allow_dt_increase) throw oxDNAException ("should never happen");

	if (fact < 1.f) fact = 1.f;
	printf ("halfway fact: %g\n", fact);

	// here we increment coordinates
	for (int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->pos += p->vel * (this->_dt / fact);
		
		if (p->is_rigid_body()) {
			LR_matrix<number> myR = dRs.back();
			dRs.pop_back();
			p->orientation = (p->orientation * myR)  / fact;
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
			p->torque = LR_vector<number>((number) 0, (number) 0, (number) 0);

		}
		p->set_initial_forces(cur_step, this->_box);

		this->_lists->single_update(p);
	}

	if(is_warning) {
		std::stringstream ss;
		for(vector<int>::iterator it = w_ps.begin(); it != w_ps.end(); it++) ss << *it << " ";
		OX_LOG(Logger::LOG_WARNING, "The following particles had a displacement greater than 0.1 in this step: %s", ss.str().c_str());
	}
}


template<typename number>
void FIREBackend<number>::_second_step() {
	this->_K = (number) 0.f;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->vel += p->force * this->_dt * (number) 0.5f;
		if(p->is_rigid_body()) p->L += p->torque * this->_dt * (number) 0.5f;

		this->_K += (p->vel.norm() + p->L.norm()) * (number) 0.5f;
	}
}

template<typename number>
void FIREBackend<number>::_compute_forces() {
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
void FIREBackend<number>::_evolve () {
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

	if (fact_r < 1.f) fact_r = 1.f;
	if (fact_l < 1.f) fact_l = 1.f;

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
void FIREBackend<number>::sim_step(llint curr_step) {
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

	// calculate P
	number P = 0.;
	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> * p = this->_particles[i];
		P += p->force * p->vel;
		if (p->is_rigid_body()) P += p->L * p->torque;
	}
	
	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> * p = this->_particles[i];
		p->vel = p->vel * (1. - _alpha) + p->force * (_alpha * p->vel.module() / p->force.module());
		if (p->is_rigid_body()) {
			p->L = p->L * (1. - _alpha) + p->torque * (_alpha * p->L.module() / p->torque.module());
		}
	}

	_step_last_reset += 1;
	
	if (P > (number)0.f && _step_last_reset > _N_min && _allow_dt_increase) {
		this->_dt = this->_dt * _f_inc;
		if (this->_dt > _dt_max) this->_dt = _dt_max;
		_alpha = _alpha * _f_alpha;
	}
	else if (P <= (number) 0.f) {
		this->_dt = this->_dt * _f_dec;
		_alpha = _alpha_start;
		for (int i = 0; i < this->_N; i ++) {
			BaseParticle<number> * p = this->_particles[i];
			p->vel.x = (number) 0.f;
			p->vel.y = (number) 0.f;
			p->vel.z = (number) 0.f;
			if (p->is_rigid_body()) {
				p->L.x = (number) 0.f;
				p->L.y = (number) 0.f;
				p->L.z = (number) 0.f;
			}
		}
		_step_last_reset = 0;
	}
	this->_timer_thermostat->pause();

	this->_mytimer->pause();
	printf("P, alpha, dt = %g, %g, %g\n", P, _alpha, this->_dt);
}

template class FIREBackend<float>;
template class FIREBackend<double>;


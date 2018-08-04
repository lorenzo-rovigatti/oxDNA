/*
 * MD_CPUBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MD_CPUBackend.h"
#include "./Thermostats/ThermostatFactory.h"
#include "MCMoves/MoveFactory.h"

template<typename number>
MD_CPUBackend<number>::MD_CPUBackend() : MDBackend<number>() {
	this->_is_CUDA_sim = false;
	_thermostat = NULL;
	_V_move = NULL;
	_compute_stress_tensor = false;
	_stress_tensor_avg_every = -1;
	_stress_tensor_counter = 0;
}

template<typename number>
MD_CPUBackend<number>::~MD_CPUBackend() {
	delete _thermostat;
	if(this->_use_barostat) delete _V_move;
}

template<typename number>
void MD_CPUBackend<number>::_first_step(llint curr_step) {
	bool is_warning = false;
	std::vector<int> w_ps;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];

		p->vel += p->force*(this->_dt*(number)0.5);
		LR_vector<number> dr = p->vel*this->_dt;
		if(dr.norm() > 0.01) {
			is_warning = true;
			w_ps.push_back(p->index);
		}
		p->pos += dr;
		// if Lees-Edwards boundaries are enabled, we have to check for crossings along the y axis
		if(this->_lees_edwards) {
			const LR_vector<number> &L = this->_box->box_sides();
			int y_new = floor(p->pos.y / L.y);
			int y_old = floor((p->pos.y - dr.y)/ L.y);
			// we crossed the boundary along y
			if(y_new != y_old) {
				number delta_x = this->_shear_rate * L.y * curr_step * this->_dt;
				delta_x -= floor(delta_x/L.x)*L.x;
				if(y_new > y_old) {
					p->pos.x -= delta_x;
					p->pos.y -= L.y;
					p->vel.x -= this->_shear_rate*L.y;
				}
				else {
					p->pos.x += delta_x;
					p->pos.y += L.y;
					p->vel.x += this->_shear_rate*L.y;
				}
			}
		}

		if(p->is_rigid_body()) {
			p->L += p->torque*(this->_dt*(number)0.5);
			// update of the orientation
			number norm = p->L.module();
			LR_vector<number> LVersor(p->L/norm);

			number sintheta = sin(this->_dt*norm);
			number costheta = cos(this->_dt*norm);
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

			p->orientation = p->orientation*R;
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
		}

		p->set_initial_forces(curr_step, this->_box);

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

		typename vector<ParticlePair<number> >::iterator it = p->affected.begin();
		for(; it != p->affected.end(); it++) {
			if(it->first == p) this->_U += this->_interaction->pair_interaction_bonded(it->first, it->second, NULL, true);
		}

		std::vector<BaseParticle<number> *> neighs = this->_lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle<number> *q = neighs[n];
			if(!_compute_stress_tensor)	this->_U += this->_interaction->pair_interaction_nonbonded(p, q, NULL, true);
			else _update_forces_and_stress_tensor(p, q);
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

		if(_compute_stress_tensor) _update_kinetic_stress_tensor(p);
	}
}

template<typename number>
void MD_CPUBackend<number>::_update_forces_and_stress_tensor(BaseParticle<number> *p, BaseParticle<number> *q) {
	// pair_interaction_nonbonded will change these vectors, but we still need them in the next
	// first integration step. For this reason we copy and then restore their values
	// after the calculation
	LR_vector<number> old_p_force(p->force);
	LR_vector<number> old_q_force(q->force);
	LR_vector<number> old_p_torque(p->torque);
	LR_vector<number> old_q_torque(q->torque);

	p->force = q->force = p->torque = q->torque = LR_vector<number>();

	LR_vector<number> r = this->_box->min_image(p->pos, q->pos);
	this->_interaction->pair_interaction_nonbonded(p, q, &r, true);

	_stress_tensor.v1.x += r.x * q->force.x;
	_stress_tensor.v1.y += r.x * q->force.y;
	_stress_tensor.v1.z += r.x * q->force.z;
	_stress_tensor.v2.x += r.y * q->force.x;
	_stress_tensor.v2.y += r.y * q->force.y;
	_stress_tensor.v2.z += r.y * q->force.z;
	_stress_tensor.v3.x += r.z * q->force.x;
	_stress_tensor.v3.y += r.z * q->force.y;
	_stress_tensor.v3.z += r.z * q->force.z;

	p->force += old_p_force;
	q->force += old_q_force;
	p->torque += old_p_torque;
	q->torque += old_q_torque;
}

template<typename number>
void MD_CPUBackend<number>::_update_kinetic_stress_tensor(BaseParticle<number> *p) {
	LR_vector<number> vel = p->vel;
	if(this->_shear_rate > 0.) {
		number Ly = this->_box->box_sides().y;
		number y_in_box = p->pos.y - floor(p->pos.y/Ly)*Ly - 0.5*Ly;
		number flow_vx = y_in_box*this->_shear_rate;
		vel.x -= flow_vx;
	}
	_stress_tensor.v1.x += SQR(vel.x);
	_stress_tensor.v1.y += vel.x * vel.y;
	_stress_tensor.v1.z += vel.x * vel.z;
	_stress_tensor.v2.x += vel.y * vel.x;
	_stress_tensor.v2.y += SQR(vel.y);
	_stress_tensor.v2.z += vel.y * vel.z;
	_stress_tensor.v3.x += vel.z * vel.x;
	_stress_tensor.v3.y += vel.z * vel.y;
	_stress_tensor.v3.z += SQR(vel.z);
}

template<typename number>
void MD_CPUBackend<number>::_update_backend_info() {
	_stress_tensor /= 3*this->_box->V();
	double P = _stress_tensor.v1.x + _stress_tensor.v2.y + _stress_tensor.v3.z;
	this->_backend_info = Utils::sformat("% 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf", P, _stress_tensor.v1.x, _stress_tensor.v1.y, _stress_tensor.v1.z, _stress_tensor.v2.x, _stress_tensor.v2.y, _stress_tensor.v2.z, _stress_tensor.v3.x, _stress_tensor.v3.y, _stress_tensor.v3.z);
	_stress_tensor = LR_matrix<double>();
}

template<typename number>
void MD_CPUBackend<number>::sim_step(llint curr_step) {
	this->_mytimer->resume();

	CONFIG_INFO->curr_step = curr_step;

	this->_timer_first_step->resume();
	_first_step(curr_step);
	this->_timer_first_step->pause();

	this->_timer_lists->resume();
	if(!this->_lists->is_updated()) {
		this->_lists->global_update();
		this->_N_updates++;
	}
	this->_timer_lists->pause();

	if(this->_is_barostat_active()) {
		this->_timer_barostat->resume();
		_V_move->apply(curr_step);
		this->_barostat_acceptance = _V_move->get_acceptance();
		this->_timer_barostat->pause();
	}

	this->_timer_forces->resume();
	_compute_forces();
	_second_step();

	if(_compute_stress_tensor) {
		_stress_tensor_counter++;
		if(_stress_tensor_counter == _stress_tensor_avg_every) {
			_stress_tensor_counter = 0;
			_stress_tensor /= _stress_tensor_avg_every;
			_update_backend_info();
		}
	}

	this->_timer_forces->pause();

	this->_timer_thermostat->resume();
	_thermostat->apply(this->_particles, curr_step);
	this->_timer_thermostat->pause();

	this->_mytimer->pause();
}

template<typename number>
void MD_CPUBackend<number>::get_settings (input_file &inp) {
	MDBackend<number>::get_settings(inp);

	_thermostat = ThermostatFactory::make_thermostat<number>(inp, this->_box);
	_thermostat->get_settings(inp);

	getInputBool(&inp, "MD_compute_stress_tensor", &_compute_stress_tensor, 0);
	if(_compute_stress_tensor) {
		OX_LOG(Logger::LOG_INFO, "Computing the stress tensor directly in the backend");
		getInputInt(&inp, "MD_stress_tensor_avg_every", &_stress_tensor_avg_every, 1);
		if(_stress_tensor_avg_every < 1) throw oxDNAException("MD_stress_tensor_avg_every should be larger than 0");
	}

	if(this->_use_barostat) {
		string str_inp = Utils::sformat("type = volume\ndelta = %lf\nisotropic = %d", this->_delta_L, (int) this->_barostat_isotropic);
		input_file *move_inp = Utils::get_input_file_from_string(str_inp);
		_V_move = MoveFactory::make_move<number>(*move_inp, inp);
		cleanInputFile(move_inp);
	}
}

template<typename number>
void MD_CPUBackend<number>::init() {
	MDBackend<number>::init();
	_thermostat->init (this->_N);
	if(this->_use_barostat) _V_move->init();

	_compute_forces();
	if(_compute_stress_tensor) {
		for(int i = 0; i < this->_N; i++) {
			BaseParticle<number> *p = this->_particles[i];
			_update_kinetic_stress_tensor(p);
		}
		_update_backend_info();
	}
}

template class MD_CPUBackend<float>;
template class MD_CPUBackend<double>;

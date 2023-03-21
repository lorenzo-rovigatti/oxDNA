/*
 * Pivot.cpp
 *
 *  Created on: Mar 8, 2023
 *      Author: lorenzo
 */

#include "Pivot.h"

Pivot::Pivot() :
				BaseMove() {

}

Pivot::~Pivot() {

}

void Pivot::get_settings(input_file &inp, input_file &sim_inp) {
	BaseMove::get_settings(inp, sim_inp);

	getInputNumber(&inp, "delta", &_delta, 1);
}

void Pivot::init() {
	BaseMove::init();
	OX_LOG(Logger::LOG_INFO, "(MCRot.cpp) Pivot initiated with T %g, delta %g, prob: %g", _T, _delta, prob);
}

void Pivot::_copy_particle_details(BaseParticle *dest, BaseParticle *src) {
	dest->index = src->index;
	dest->pos = src->pos;
	dest->orientation = src->orientation;
	dest->orientationT = src->orientationT;
}

void Pivot::apply(llint curr_step) {
	_attempted++;

	static std::vector<BaseParticle> temp_particles(_Info->N());
	if(temp_particles.size() > (uint) _Info->N()) {
		temp_particles.resize(_Info->N());
	}

	int pivot_idx = (int) (drand48() * _Info->N());
	BaseParticle *pivot_p = _Info->particles()[pivot_idx];

	// pick a direction and follow it till we reach the end of the chain
	int dir = (drand48() > 0.5) ? +1 : -1;
	BaseParticle *next_p = pivot_p;
	int N_in_move = 0;
	bool done = false;
	do {
		int next_idx = next_p->index + dir;
		if(next_idx >= 0 && next_idx < _Info->N()) {
			next_p = _Info->particles()[next_idx];
			if(next_p->strand_id == pivot_p->strand_id) {
				_copy_particle_details(&temp_particles[N_in_move], next_p);
				N_in_move++;
			}
			else {
				done = true;
			}
		}
		else {
			done = true;
		}
	}
	while(!done);

	// the section of the chain that should be moved is empty -> early rejection
	if(N_in_move == 0) {
		return;
	}

	// compute the energy in the initial state
	number delta_E = 0.;
	for(int i = 0; i < N_in_move; i++) {
		BaseParticle *p = _Info->particles()[temp_particles[i].index];
		p->set_ext_potential(curr_step, _Info->box);
		delta_E -= particle_energy(p) + p->ext_potential;
	}

	// move all the particles involved
	LR_matrix R = Utils::get_random_rotation_matrix_from_angle(drand48() * _delta);
	for(int i = 0; i < N_in_move; i++) {
		BaseParticle *p = _Info->particles()[temp_particles[i].index];
		p->pos = R * (p->pos - pivot_p->pos) + pivot_p->pos;
		p->orientation = p->orientation * R;
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();
		_Info->lists->single_update(p);
	}

	if(!_Info->lists->is_updated()) {
		_Info->lists->global_update();
	}

	// compute the energy in the final state
	for(int i = 0; i < N_in_move; i++) {
		BaseParticle *p = _Info->particles()[temp_particles[i].index];
		p->set_ext_potential(curr_step, this->_Info->box);
		delta_E += particle_energy(p) + p->ext_potential;
	}

	// accept or reject?
	if(_Info->interaction->get_is_infinite() == false && (delta_E  < 0 || exp(-delta_E / _T) > drand48() )) {
		// move accepted
		_accepted++;
		if (curr_step < _equilibration_steps && _adjust_moves) {
			_delta *= _acc_fact;
			if(_delta > M_PI) {
				_delta = M_PI;
			}
		}

		_Info->lists->global_update();
	}
	else {
		for(int i = 0; i < N_in_move; i++) {
			BaseParticle *p = _Info->particles()[temp_particles[i].index];
			_copy_particle_details(p, &temp_particles[i]);
			p->set_positions();
			_Info->lists->single_update(p);
		}

		if(!_Info->lists->is_updated()) {
			_Info->lists->global_update();
		}

		for(int i = 0; i < N_in_move; i++) {
			BaseParticle *p = _Info->particles()[temp_particles[i].index];
			p->set_ext_potential(curr_step, _Info->box);
		}

		_Info->interaction->set_is_infinite(false);

		if(curr_step < _equilibration_steps && _adjust_moves) {
			_delta /= _rej_fact;
		}
	}
}

void Pivot::log_parameters() {
	BaseMove::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}

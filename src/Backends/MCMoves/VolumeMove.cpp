/**
 * @file    VolumeMove.cpp
 * @date    30/may/2015
 * @author  flavio
 *
 *
 */

#include "VolumeMove.h"

/// traslation
template<typename number>
VolumeMove<number>::VolumeMove (ConfigInfo<number> * Info) : BaseMove<number>(Info) {
	_verlet_skin = -1.f;
}

template<typename number>
VolumeMove<number>::~VolumeMove () {

}

template<typename number>
void VolumeMove<number>::init () {
	BaseMove<number>::init();
	_pos_old.resize (*this->_Info->N);
}

template<typename number>
void VolumeMove<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "prob", &this->prob, 0);
	getInputNumber (&sim_inp, "P", &_P, 1);
	OX_LOG(Logger::LOG_INFO, "(VolumeMove.cpp) VolumeMove initiated with T %g, delta %g, prob: %g", this->_T, _delta, this->prob);
	
	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}

template<typename number>
void VolumeMove<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	BaseParticle<number> ** particles = this->_Info->particles;
	int N = *(this->_Info->N);
	
	number dL = _delta * (drand48() - (number) 0.5);
	LR_vector<number> box_sides = this->_Info->box->box_sides();
	LR_vector<number> old_box_sides = this->_Info->box->box_sides();

	number oldE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	number oldV = (box_sides[0]) * (box_sides[0]) * (box_sides[0]);

	this->_Info->interaction->set_box_side(box_sides[0] + dL);
	this->_Info->box->init(box_sides[0] + dL, box_sides[1] + dL, box_sides[2] + dL);
	this->_Info->lists->change_box();
	number dExt = (number) 0.f;
	for (int k = 0; k < N; k ++) {
		BaseParticle<number> *p = particles[k];
		dExt = -p->ext_potential;
		_pos_old[k] = p->pos; 
		p->pos *= (1. + dL / box_sides[0]);
		p->set_ext_potential(curr_step, box_sides[0]);
		dExt += -p->ext_potential;
	}

	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	number newE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	number dE = newE - oldE + dExt;
	number V = (box_sides[0] + dL) * (box_sides[0] + dL) * (box_sides[0] + dL);
	number dV = V - oldV;

	if (this->_Info->interaction->get_is_infinite() == false &&
			exp (-(dE + _P * dV - N * this->_T * log (V / oldV)) / this->_T) > drand48()) {
		this->_accepted ++;
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta *= this->_acc_fact;
		*(this->_Info->box_side) = box_sides[0];
	}
	else {
		for (int k = 0; k < N; k ++) {
			BaseParticle<number> *p = particles[k];
			p->pos = _pos_old[k]; 
			p->set_ext_potential(curr_step, old_box_sides[0]);
		}
		this->_Info->interaction->set_box_side(old_box_sides[0]);
		this->_Info->box->init(old_box_sides[0], old_box_sides[1], old_box_sides[2]);
		this->_Info->lists->change_box();
		this->_Info->interaction->set_is_infinite (false);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}
	return;
}

template class VolumeMove<float>;
template class VolumeMove<double>;

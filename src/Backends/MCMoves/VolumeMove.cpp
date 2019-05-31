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
VolumeMove<number>::VolumeMove ()  {
	_verlet_skin = -1.f;
	_isotropic = true;
}

template<typename number>
VolumeMove<number>::~VolumeMove () {

}

template<typename number>
void VolumeMove<number>::init () {
	BaseMove<number>::init();
	_pos_old.resize (*this->_Info->N);
	if (this->_restrict_to_type > 0) OX_LOG (Logger::LOG_WARNING, "(VolumeMove.cpp) Cant use VolumeMove with restrict_to_type. Ignoring");
	OX_LOG(Logger::LOG_INFO, "(VolumeMove.cpp) VolumeMove (isotropic = %d) initiated with T %g, delta %g, prob: %g", _isotropic, this->_T, _delta, this->prob);
}

template<typename number>
void VolumeMove<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputBool(&inp, "isotropic", &_isotropic, 0);
	getInputNumber(&inp, "delta", &_delta, 1);
	getInputNumber(&inp, "prob", &this->prob, 0);
	getInputNumber(&sim_inp, "P", &_P, 1);

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

	LR_vector<number> box_sides = this->_Info->box->box_sides();
	LR_vector<number> old_box_sides = box_sides;

	number oldE;
	if (this->_compute_energy_before) oldE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	else oldE = (number) 0.f;
	number oldV = this->_Info->box->V();

	if(_isotropic) {
		number dL = _delta*(drand48() - (number) 0.5);
		box_sides.x += dL;
		box_sides.y += dL;
		box_sides.z += dL;
		/*number newL = pow(exp(log(oldV) + _delta*(drand48() - (number) 0.5)), 1./3.);
		  box_sides.x = box_sides.y = box_sides.z = newL;*/
	}
	else {
		box_sides.x += _delta*(drand48() - (number) 0.5);
		box_sides.y += _delta*(drand48() - (number) 0.5);
		box_sides.z += _delta*(drand48() - (number) 0.5);
	}

	this->_Info->box->init(box_sides[0], box_sides[1], box_sides[2]);
	number dExt = (number) 0.f;
	for(int k = 0; k < N; k ++) {
		BaseParticle<number> *p = particles[k];
		dExt -= p->ext_potential;
		_pos_old[k] = p->pos;
		p->pos.x *= box_sides[0]/old_box_sides[0];
		p->pos.y *= box_sides[1]/old_box_sides[1];
		p->pos.z *= box_sides[2]/old_box_sides[2];
		p->set_ext_potential(curr_step, this->_Info->box);
		dExt += p->ext_potential;
	}

	// this bit has to come after the update of particles' positions
	this->_Info->lists->change_box();
	if(!this->_Info->lists->is_updated()) this->_Info->lists->global_update();

	number newE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	number dE = newE - oldE + dExt;
	number V = this->_Info->box->V();
	number dV = V - oldV;

	if (this->_Info->interaction->get_is_infinite() == false && exp(-(dE + _P*dV - (N + 1)*this->_T*log(V/oldV))/this->_T) > drand48()) {
		//printf("B %lld %lf ---- %lf %lf %lf\n", curr_step, newE/N, dE, dV, exp(-(dE + _P*dV - N*this->_T*log(V/oldV))/this->_T));
		this->_accepted++;
		if(curr_step < this->_equilibration_steps && this->_adjust_moves) _delta *= this->_acc_fact;
	}
	else {
		for (int k = 0; k < N; k ++) {
			BaseParticle<number> *p = particles[k];
			p->pos = _pos_old[k];
			p->set_ext_potential(curr_step, this->_Info->box);
		}
		this->_Info->interaction->set_is_infinite(false);
		this->_Info->box->init(old_box_sides[0], old_box_sides[1], old_box_sides[2]);
		this->_Info->lists->change_box();
		if(!this->_Info->lists->is_updated()) this->_Info->lists->global_update();
		if(curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}
	return;
}

template<typename number>
void VolumeMove<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta", _delta);
}

template class VolumeMove<float>;
template class VolumeMove<double>;

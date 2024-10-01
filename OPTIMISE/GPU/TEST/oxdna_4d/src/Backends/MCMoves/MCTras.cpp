/**
 * @file    MCTras.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#include "MCTras.h"

/// traslation

MCTras::MCTras()  {
	pos_old = LR_vector (0., 0., 0.);

	_verlet_skin = -1.f;
}


MCTras::~MCTras () {

}


void MCTras::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);

	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}


void MCTras::init() {
	BaseMove::init();
	OX_LOG(Logger::LOG_INFO, "(MCTras.cpp) MCtras initiated with T %g, delta %g, prob: %g", this->_T, _delta, this->prob);
}


void MCTras::apply(llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * this->_Info->N());
	BaseParticle *p = this->_Info->particles()[pi];
	if (this->_restrict_to_type >= 0) {
		while(p->type != this->_restrict_to_type) {
			pi = (int) (drand48() * this->_Info->N());
			p = this->_Info->particles()[pi];
		}
	}

	pos_old = p->pos;

	// compute the energy before the move
	number delta_E;
	if (this->_compute_energy_before) delta_E = -this->particle_energy(p);
	else delta_E = (number) 0.f;
	p->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential;

	// perform the move
	p->pos.x += 2. * (drand48() - (number)0.5f) * _delta;
	p->pos.y += 2. * (drand48() - (number)0.5f) * _delta;
	p->pos.z += 2. * (drand48() - (number)0.5f) * _delta;

	// update lists
	this->_Info->lists->single_update(p);
	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	// energy after the move
	delta_E += this->particle_energy(p);
	p->set_ext_potential(curr_step, this->_Info->box);
	delta_E_ext += p->ext_potential;

	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48() )) {
		// move accepted
		this->_accepted ++;
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta *= this->_acc_fact;
			if (_delta > 0.5) _delta = 0.5;
			if (_delta < 1.e-5) _delta = 1.e-5;
			if (_verlet_skin > 0. && _delta > _verlet_skin * 0.8660254037844386 / 2.) {
				_delta = _verlet_skin *_verlet_skin * 0.8660254037844386 / 2.; //  0.8660254 .. ==sqrt(3.) / 2.
			}
		}
	}
	else {
		this->_Info->particles()[pi]->pos = pos_old;
		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);

		if(curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}

	return;
}


void MCTras::log_parameters() {
	BaseMove::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}

/**
 * @file    Reappear.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 */
#include "Reappear.h"

/// traslation
Reappear::Reappear (){
	pos_old = LR_vector (0., 0., 0.);
}

Reappear::~Reappear () {

}

void Reappear::init () {
	BaseMove::init();
	OX_LOG(Logger::LOG_INFO, "(Reappear.cpp) Reappear move initiated with probability %g", this->prob);
}

void Reappear::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove::get_settings (inp, sim_inp);
}


void Reappear::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * (_Info->N()));
	BaseParticle * p = _Info->particles()[pi];
	if (this->_restrict_to_type >= 0) {
		while(p->type != this->_restrict_to_type) {
			pi = (int) (drand48() * (_Info->N()));
			p = _Info->particles()[pi];
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
	LR_vector sides = this->_Info->box->box_sides();
	p->pos.x = drand48() * sides[0];
	p->pos.y = drand48() * sides[1];
	p->pos.z = drand48() * sides[2];

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

	}
	else {
		_Info->particles()[pi]->pos = pos_old;
		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);
	}

	return;
}


/**
 * @file    Swim.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 */
#include "Swim.h"

/// traslation
template<typename number>
Swim<number>::Swim (){
	pos_old = LR_vector<number> (0., 0., 0.);
	_delta = (number) -1.f;
	_axis_i = -1;
}

template<typename number>
Swim<number>::~Swim () {

}

template<typename number>
void Swim<number>::init () {
	BaseMove<number>::init();
	OX_LOG(Logger::LOG_INFO, "(Swim.cpp) Swim move initiated on axis %d with delta %g and probability %g", _axis_i, _delta, this->prob);
}

template<typename number>
void Swim<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputInt (&inp, "axis", &_axis_i, 1);
}

template<typename number>
LR_vector<number> * Swim<number>::_get_axis(BaseParticle<number> * p) {
	switch (_axis_i) {
		case Swim::O_V1: return &(p->orientation.v1);
		case Swim::O_V2: return &(p->orientation.v2);
		case Swim::O_V3: return &(p->orientation.v3);
		default: throw oxDNAException ("Should never get here....");
	}
	return (LR_vector<number> *) NULL;
}

template<typename number>
void Swim<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * p = this->_Info->particles[pi];
	if (this->_restrict_to_type >= 0) {
		while(p->type != this->_restrict_to_type) {
			pi = (int) (drand48() * (*this->_Info->N));
			p = this->_Info->particles[pi];
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
	p->pos += (number)((drand48() - (number)0.5f) * _delta) * (*_get_axis(p));

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

		if (curr_step < this->_equilibration_steps && this->_adjust_moves)
			_delta *= this->_acc_fact;
	}
	else {
		this->_Info->particles[pi]->pos = pos_old;

		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);

		if (curr_step < this->_equilibration_steps && this->_adjust_moves)
			_delta /= this->_rej_fact;
	}

	return;
}

template<typename number>
void Swim<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}

template class Swim<float>;
template class Swim<double>;

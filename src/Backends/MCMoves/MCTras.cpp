/**
 * @file    MCTras.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#include "MCTras.h"

/// traslation
template<typename number>
MCTras<number>::MCTras (ConfigInfo<number> * Info) : BaseMove<number>(Info) {
	pos_old = LR_vector<number> (0., 0., 0.);

	_verlet_skin = -1.f;
}

template<typename number>
MCTras<number>::~MCTras () {

}

template<typename number>
void MCTras<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	OX_LOG(Logger::LOG_INFO, "(MCTras.cpp) MCtras initiated with T %g, delta %g, prob: %g", this->_T, _delta, this->prob);
	
	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}

template<typename number>
void MCTras<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * (*this->_Info->N));

	//cout << "GGB " << this->_Info->particles << endl;;

	BaseParticle<number> *p = this->_Info->particles[pi];

	pos_old = p->pos; 

	// compute the energy before the move
	number delta_E = -this->particle_energy(p);
	p->set_ext_potential (curr_step, *this->_Info->box_side);
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
	p->set_ext_potential(curr_step, *this->_Info->box_side);
	delta_E_ext += p->ext_potential;

	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48() )) {
		// move accepted
		this->_accepted ++;
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta *= this->_acc_fact;
			if (_verlet_skin > 0. && _delta > _verlet_skin * 0.8660254037844386 / 2.) {
				_delta = _verlet_skin *_verlet_skin * 0.8660254037844386 / 2.; //  0.8660254 .. ==sqrt(3.) / 2.
			}
		}
	}
	else {
		this->_Info->particles[pi]->pos = pos_old;
		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}

	return;
}

template class MCTras<float>;
template class MCTras<double>;

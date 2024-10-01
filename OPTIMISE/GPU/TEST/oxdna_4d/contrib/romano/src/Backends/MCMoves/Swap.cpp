/**
 * @file    Swap.cpp
 * @date    10/lug/2019
 * @author  flavio
 *
 */
#include "Swap.h"

template<typename number>
Swap<number>::Swap (){

}

template<typename number>
Swap<number>::~Swap () {

}

template<typename number>
void Swap<number>::init () {
	BaseMove<number>::init();
	OX_LOG(Logger::LOG_INFO, "(Swap.cpp) Swap move initiated with probability %g", this->prob);
}

template<typename number>
void Swap<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);
}

template<typename number>
inline void do_swap(BaseParticle<number> *p, BaseParticle<number> *q) {
	// exchange pos
	LR_vector<number> tmp_pos = p->pos;
	p->pos = q->pos;
	q->pos = tmp_pos;

	// importantly, we don't exchange type: we tranlate red particle
	// into position and orientation of blue particle and vice versa

	// exchange orientation matrixes
	LR_matrix<number> tmpm = p->orientation;
	p->orientation = q->orientation;
	q->orientation = tmpm;
	tmpm = p->orientationT;
	p->orientationT = q->orientationT;
	q->orientationT = tmpm;
	
	// set the positions
	p->set_positions();
	q->set_positions();
}

template<typename number>
void Swap<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * (*this->_Info->N));
	int qi = pi;
	while (qi == pi) qi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * p = this->_Info->particles[pi];
	BaseParticle<number> * q = this->_Info->particles[qi];

	// FIXME: does it make sense to restrict to different types?
	if (p->type == q->type)
		return;

	//pos_old_p = p->pos;
	//pos_old_q = q->pos;

	// compute the energy before the move
	number delta_E = (number) 0.f;
	if (this->_compute_energy_before) {
		delta_E -= this->particle_energy(p);
		delta_E -= this->particle_energy(q);
		delta_E += this->_Info->interaction->pair_interaction(p, q); // if p and q do interact, this interaction is counted twice in the above rows
	}
	else {
		delta_E += this->_Info->interaction->pair_interaction(p, q);
	}
	p->set_ext_potential (curr_step, this->_Info->box);
	q->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential - q->ext_potential;

	// perform the move
	do_swap(p, q);

	// update lists
	this->_Info->lists->single_update(p);
	this->_Info->lists->single_update(q);
	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	// energy after the move
	delta_E += this->particle_energy(p);
	delta_E += this->particle_energy(q);
	delta_E -= this->_Info->interaction->pair_interaction(p, q);
	p->set_ext_potential(curr_step, this->_Info->box);
	q->set_ext_potential(curr_step, this->_Info->box);
	delta_E_ext += p->ext_potential + q->ext_potential;

	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48() )) {
		// move accepted
        this->_accepted ++;
	}
	else {
		do_swap(q, p);
		this->_Info->lists->single_update(p);
		this->_Info->lists->single_update(q);
		this->_Info->interaction->set_is_infinite(false);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
	}

	return;
}

template class Swap<float>;
template class Swap<double>;

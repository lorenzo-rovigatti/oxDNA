/**
 * @file    Grow.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 */
#include "Grow.h"

/// traslation
template<typename number>
Grow<number>::Grow (){
	_target_sigma = 0.5f;
	_n_to_grow = 0;
}

template<typename number>
Grow<number>::~Grow () {

}

template<typename number>
void Grow<number>::init () {
	BaseMove<number>::init();
	_particles_to_grow = std::set<int> ();
	for (int k = 0; k < *this->_Info->N; k++) {
		BaseParticle<number> * p = this->_Info->particles[k];
		if (this->_restrict_to_type >= 0) {
			if (p->type == this->_restrict_to_type) {
				_particles_to_grow.insert(k);
			}
		}
		else {
			_particles_to_grow.insert(k);
		}
	}
	_n_to_grow = _particles_to_grow.size();
	OX_LOG(Logger::LOG_INFO, "(Grow.cpp) Grow move initiated with probability %g and target_sigma=%g", this->prob, _target_sigma);
	OX_LOG(Logger::LOG_INFO, "(Grow.cpp)                will have to grow %d particles", _n_to_grow);
}

template<typename number>
void Grow<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);
	getInputNumber(&inp, "target_sigma", &_target_sigma, 1);
}


template<typename number>
void Grow<number>::apply (llint curr_step) {

	if (_particles_to_grow.size() == 0) return;

	//if (_n_grown == *this->_Info->N) return;

	//if (_n_grown == *this->_Info->N) throw oxDNAException ("DONE!");

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to grow
	int pi_h = (int) (lrand48() % _particles_to_grow.size());

	std::set<int>::iterator pit = _particles_to_grow.begin();
	//it = it + (size_t)pi_h;
	std::advance(pit, pi_h);
	if (pit == _particles_to_grow.end()) throw oxDNAException("merda");
	int pi = *pit;
	/*
	int pi = (int) (drand48() * (_particles_to_grow.size()));
	if (this->_restrict_to_type >= 0) {
		while(this->_Info->particles[pi]->type != this->_restrict_to_type || _grown[pi] == true) {
			pi = (int) (drand48() * (*this->_Info->N));
		}
	}
	*/
	SpheroCylinder<number> * s = static_cast< SpheroCylinder<number> *>(this->_Info->particles[pi]);
	BaseParticle<number> * p = this->_Info->particles[pi];

	number old_length = s->length;

	// find the minimum distance
	std::vector<BaseParticle<number> *> neighs = this->_Info->lists->get_complete_neigh_list(p);
	typename std::vector<BaseParticle<number> *>::iterator it;
	number rmin = _target_sigma;
	for (it = neighs.begin(); it != neighs.end(); it ++) {
		number myr = 1. + rmin;
		SpheroCylinder<number> * q;
		if (this->_restrict_to_type >= 0) {
			if ((*it)->type == this->_restrict_to_type) {
				q = static_cast<SpheroCylinder<number> *> (*it);
				if (p->index < q->index) myr = sqrt(this->_Info->box->sqr_min_image_distance(p, *it)) - q->length;
				else myr = sqrt(this->_Info->box->sqr_min_image_distance(*it, p)) - q->length;
			}
		}
		else {
			q = static_cast<SpheroCylinder<number> *> (*it);
			if (p->index < q->index) myr = sqrt(this->_Info->box->sqr_min_image_distance(p, *it)) - q->length;
			else myr = sqrt(this->_Info->box->sqr_min_image_distance(*it, p)) - q->length;
		}
		if (myr < rmin) rmin = myr;
	}

	if (rmin < old_length) rmin = old_length;

	// compute the energy before the move
	number delta_E;
	if (this->_compute_energy_before) delta_E = -this->particle_energy(p);
	else delta_E = (number) 0.f;
	p->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential;

	if (rmin >= _target_sigma) s->length = _target_sigma;
	else s->length = rmin;

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
		if (s->length >= _target_sigma) _particles_to_grow.erase(pi);
		if (_particles_to_grow.size() == 0) {
			OX_LOG(Logger::LOG_INFO, "(Grow.cpp) Growing complete!");
		}
	}
	else {
		s->length = old_length;
		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);

		this->particle_energy(p);
		if (this->_Info->interaction->get_is_infinite() == true) throw oxDNAException("disaster");
	}

	if (this->_attempted % 100000 == 0) OX_LOG(Logger::LOG_INFO, "(Grow.cpp) have grown %d particles out of %d, (%5.3lf)", _n_to_grow - _particles_to_grow.size(), _n_to_grow, 1. - ((double)_particles_to_grow.size() / (double)_n_to_grow));

	return;
}

template<typename number>
double Grow<number>::get_acceptance() {
	if (true || _particles_to_grow.size() == 0) return ((double)this->_accepted)/((double) this->_attempted);
	else return 1. - ((double)_particles_to_grow.size() / (double)_n_to_grow);
}

template class Grow<float>;
template class Grow<double>;

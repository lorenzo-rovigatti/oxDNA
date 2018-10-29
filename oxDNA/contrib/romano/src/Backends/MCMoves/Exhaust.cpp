/**
 * @file    Exhaust.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 */
#include "Exhaust.h"

/// traslation
template<typename number>
Exhaust<number>::Exhaust (){
	_delta = (number) -1.f;
	_poss_old = std::vector<LR_vector<number> > ();
	_clust = std::vector<BaseParticle<number> *> ();
	_max_clust_size = -1;
	_n_found_other = (llint) 0;
	_n_too_large = (llint) 0;
}

template<typename number>
Exhaust<number>::~Exhaust () {

}

template<typename number>
void Exhaust<number>::init () {
	BaseMove<number>::init();
	if (this->_restrict_to_type < 0) throw oxDNAException("Exhaust move MUST be restricted to a type");
	_clust.reserve((size_t) _max_clust_size);
	_poss_old.reserve((size_t) _max_clust_size);
	OX_LOG(Logger::LOG_INFO, "(Exhaust.cpp) Exhaust move initiated with delta %g, max_clust_size=%d and probability %g", _delta, _max_clust_size, this->prob);
}

template<typename number>
void Exhaust<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);
	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "delta_max", &_delta_max, 1);
	getInputInt (&inp, "max_clust_size", &_max_clust_size, 1);
}

template<typename number>
void Exhaust<number>::apply (llint curr_step) {

	_clust.resize(0);
	_poss_old.resize(0);

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

	// compute the energy before the move
	number delta_E;
	if (this->_compute_energy_before) delta_E = -this->particle_energy(p);
	else delta_E = (number) 0.f;
	p->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential;

	LR_vector<number> disp_dir = Utils::get_random_vector<number>();
	//number disp_mag = drand48() * _delta;

	typename std::vector<BaseParticle<number> *>::iterator it;

	// find closest particle to bump into
	_clust.push_back(p);
	_poss_old.push_back(p->pos);
	bool recruited = true;
	bool found_other = false;
	number remaining = _delta;
	while (recruited == true && (int)_clust.size() <= _max_clust_size) {
		recruited = false;
		std::vector<BaseParticle<number> *> neighs = this->_Info->lists->get_complete_neigh_list(p);

		number lambdamin = _delta;
		int idxmin = -1;

		p = _clust.back();
		//printf ("starting from particle %d\n", p->index);

		for (it = neighs.begin(); it != neighs.end(); it ++) {
			BaseParticle<number> * q = *it;
			if (this->_restrict_to_type >= 0) {
				if ((*it)->type == this->_restrict_to_type) {
					LR_vector<number> myr = this->_Info->box->min_image(p, q);
					number proj = myr * disp_dir; // can be negative
					if (proj > (number) 0.f) {
						number myrnorm2 = myr.norm();
						number lambda;
						number det = proj*proj - (myrnorm2 - (number)1.f);
						if (det > 0.f) {
							// if det <= 0, we skip the computation
							lambda = proj - sqrt(det);

							// lambda is also the module of the displacement
							if (lambda < lambdamin) {
								lambdamin = lambda;
								idxmin = q->index;
							}
						}
					}
				}
				else {
					found_other = true;
				}
	 		}
		}

		//printf ("   lambdamin, idxmin, %g, %d ; ", lambdamin, idxmin);
		if (lambdamin < remaining && idxmin != -1) {
			BaseParticle<number> * q = this->_Info->particles[idxmin];
			_clust.push_back(q);
			_poss_old.push_back(q->pos);
			recruited = true;
			remaining = remaining - lambdamin;
		//	printf ("   recruited %d\n", _clust.back()->index);
		}
		else {
			lambdamin = remaining;
			remaining = (number) 0.f;
		//	printf ("\n");
		}
		//printf ("   moving %d by %g (%g remainig)\n", p->index, lambdamin, remaining);

		p->pos = p->pos + lambdamin * disp_dir;
		// update lists
		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
	}

	//printf ("final cluster: ");
	//for (size_t k = 0; k < _clust.size(); k++) printf(" %d", _clust[k]->index);
	//printf ("\n\n");

	//if (_clust.size() > (size_t) 4) abort();
	if (found_other) _n_found_other ++;

	//if ((int)_clust.size() <= _max_clust_size && _clust.size() != _poss_old.size()) throw oxDNAException("Should not happen %s %d (%d %d)", __FILE__, __LINE__, (int)_clust.size(), (int)_poss_old.size());
	if (_clust.size() != _poss_old.size()) throw oxDNAException("Should not happen %s %d", __FILE__, __LINE__);

	bool outright_reject = false;
	//if ((int)_clust.size() <= _max_clust_size && found_other == false) {
	if ((int) _clust.size() <= _max_clust_size) {
		for (it = _clust.begin(); it != _clust.end(); ++it) {
			p = *it;
			// energy after the move
			delta_E += this->particle_energy(p);
			p->set_ext_potential(curr_step, this->_Info->box);
			delta_E_ext += p->ext_potential;
		}
	}
	else {
		_n_too_large ++;
		outright_reject = true;
	}

	// accept or reject?
	if (outright_reject == false &&
		this->_Info->interaction->get_is_infinite() == false &&
		((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48()) ) {
		// move accepted
        this->_accepted ++;

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta *= this->_acc_fact;
			if (_delta > _delta_max) _delta = _delta_max;
		}
	}
	else {
		while (_clust.size() > 0) {
			p = _clust.back();
			p->pos = _poss_old.back();
			this->_Info->lists->single_update(p);
			this->_Info->interaction->set_is_infinite(false);
			if(!this->_Info->lists->is_updated()) {
				this->_Info->lists->global_update();
			}
			_clust.pop_back();
			_poss_old.pop_back();
		}
		if (curr_step < this->_equilibration_steps && this->_adjust_moves)
			_delta /= this->_rej_fact;
	}

	return;
}

template<typename number>
void Exhaust<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g, delta_max = %g, max_clust_size=%d (%g,%g)", _delta, _delta_max, _max_clust_size, _n_found_other / (double)this->_attempted, _n_too_large/(double) this->_attempted);
}

template class Exhaust<float>;
template class Exhaust<double>;

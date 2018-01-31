/*
 * LevyDelta.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "LevyDelta.h"

#include <sstream>

using namespace std;

template<typename number>
LevyDelta<number>::LevyDelta() {
	_tries = 0;
	_patchy_rcut = 0.16;
	_inter = NULL;
}

template<typename number>
LevyDelta<number>::~LevyDelta() {

}

template<typename number>
void LevyDelta<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "tries", &_tries, 1);
	getInputNumber(&sim_inp, "T", &_temperature, 1);
}

template<typename number>
void LevyDelta<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	_inter = dynamic_cast<LevyInteraction<number> *>(config_info.interaction);
	if(_inter == NULL) throw oxDNAException("LevyDelta can be used with LevyInteraction simulations only");

	for(int i = 0; i < *config_info.N; i++) {
		BaseParticle<number> *p = config_info.particles[i];
		if(p->btype == _inter->TETRA_CENTRE || p->btype == _inter->TETRA_PATCHY) _tetramer.push_back(p);
		else _dimer.push_back(p);
	}
}

template<typename number>
BaseParticle<number> *LevyDelta<number>::_get_random_particle(vector<BaseParticle<number> *> &v, int type) {
	BaseParticle<number> *ret;
	do {
		ret = v[lrand48() % v.size()];
	}
	while(ret->btype != type);

	return ret;
}

template<typename number>
void LevyDelta<number>::_rototranslate(vector<BaseParticle<number> *> &v , LR_vector<number> &disp, BaseParticle<number> *p, LR_matrix<number> &R) {
	for(typename vector<BaseParticle<number> *>::iterator it = v.begin(); it != v.end(); it++) {
		(*it)->pos += disp;
	}

	LR_vector<number> origin = p->pos + p->int_centers[0];
	for(typename vector<BaseParticle<number> *>::iterator it = v.begin(); it != v.end(); it++) {
		(*it)->pos = R*((*it)->pos - origin) + origin;
	}

	p->orientationT.v1 = origin - p->pos;
	Utils::orthonormalize_matrix<number>(p->orientationT);
	p->orientation = p->orientationT.get_transpose();
	p->set_positions();
}

template<typename number>
string LevyDelta<number>::get_output_string(llint curr_step) {
	stringstream outstr;

	BaseParticle<number> *rec = _get_random_particle(_tetramer, _inter->TETRA_PATCHY);
	BaseParticle<number> *p = _get_random_particle(_dimer, _inter->DIMER_PATCHY);

	number intra_energy = _inter->get_system_energy_term(LevyInteraction<number>::NONBONDED, this->_config_info.particles, *this->_config_info.N, this->_config_info.lists);

	number integrand = 0.;
	for(int i = 0; i < _tries; i++) {
		LR_vector<number> disp = rec->pos + rec->int_centers[0] - (p->pos + p->int_centers[0]) + Utils::get_random_vector_in_sphere<number>(_patchy_rcut);
		// early rejection
		LR_vector<number> dist = p->pos + disp - rec->pos;
		if(dist.norm() < 1.) continue;

		LR_matrix<number> R = Utils::get_random_rotation_matrix_from_angle<number>(2*M_PI);
		_rototranslate(_dimer, disp, p, R);
		this->_config_info.lists->global_update(true);

//		LR_vector<number> diff = p->pos + p->int_centers[0] - (rec->pos + rec->int_centers[0]);
//		printf("%lf %lf\n", diff.module(), disp.module());

		number energy = _inter->get_system_energy_term(LevyInteraction<number>::NONBONDED, this->_config_info.particles, *this->_config_info.N, this->_config_info.lists);
		number delta_energy = energy - intra_energy;
		if(delta_energy < 100) {
//			printf("%lf %lf\n", intra_energy, energy);
			integrand += exp(-delta_energy/_temperature) - 1.;
		}
	}

	LR_vector<number> disp = this->_config_info.box->box_sides()/2 + rec->pos;
	for(typename vector<BaseParticle<number> *>::iterator it = _dimer.begin(); it != _dimer.end(); it++) {
		(*it)->pos += disp;
	}
	this->_config_info.lists->global_update(true);

	integrand *= 4*M_PI*CUB(_patchy_rcut)/3.;
	outstr << integrand / (number)_tries;

	return outstr.str();
}

template class LevyDelta<float>;
template class LevyDelta<double>;

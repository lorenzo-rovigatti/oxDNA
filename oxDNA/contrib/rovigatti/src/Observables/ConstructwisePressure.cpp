/*
 * ConstructwisePressure.cpp
 *
 *  Created on: 08/jan/2018
 *      Author: lorenzo
 */

#include "ConstructwisePressure.h"

template<typename number>
ConstructwisePressure<number>::ConstructwisePressure() :
				BaseObservable<number>(),
				_P(0.),
				_shear_rate(0.),
				_construct_size(-1),
				_N_constructs(-1) {

}

template<typename number>
ConstructwisePressure<number>::~ConstructwisePressure() {

}

template<typename number>
void ConstructwisePressure<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	char raw_T[256];
	getInputString(&sim_inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number>(raw_T);

	getInputInt(&my_inp, "construct_size", &_construct_size, 1);

	bool lees_edwards = false;
	getInputBool(&sim_inp, "lees_edwards", &lees_edwards, 0);
	if(lees_edwards) getInputNumber(&sim_inp, "lees_edwards_shear_rate", &_shear_rate, 0);
}

template<typename number>
void ConstructwisePressure<number>::init(ConfigInfo<number> &info) {
	BaseObservable<number>::init(info);

	int N = *info.N;

	if(N % _construct_size) throw oxDNAException("ConstructwisePressure: the total number of particles (%d) is not a multiple of the construct size specified in the input file (%d)", N, _construct_size);
	_N_constructs = N / _construct_size;
	_construct_coms.resize(_N_constructs);
}

template<typename number>
void ConstructwisePressure<number>::update_pressure() {
	fill(_construct_coms.begin(), _construct_coms.end(), LR_vector<number>());
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		int p_construct = p->index / _construct_size;
		_construct_coms[p_construct] += this->_config_info.box->get_abs_pos(p);
	}
	for(int i = 0; i < _N_constructs; i++) {
		_construct_coms[i] /= _construct_size;
	}

	vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();

	double virial = 0;
	double energy = 0.;
	// we loop over all the pairs in order to update the forces
	typename vector<ParticlePair<number> >::iterator it;
	for(it = pairs.begin(); it != pairs.end(); it++) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;

		int p_construct = p->index / _construct_size;
		int q_construct = q->index / _construct_size;

		if(p_construct != q_construct) {
			LR_vector<number> r = this->_config_info.box->min_image(_construct_coms[p_construct], _construct_coms[q_construct]);

			// pair_interaction will change these vectors, but we still need them in the next
			// first integration step. For this reason we copy and then restore their values
			// after the calculation
			LR_vector<number> old_p_force(p->force);
			LR_vector<number> old_q_force(q->force);
			LR_vector<number> old_p_torque(p->torque);
			LR_vector<number> old_q_torque(q->torque);

			p->force = q->force = p->torque = q->torque = LR_vector<number>();

			energy += (double) this->_config_info.interaction->pair_interaction(p, q, NULL, true);

			virial -= (r * p->force);

			p->force = old_p_force;
			q->force = old_q_force;
			p->torque = old_p_torque;
			q->torque = old_q_torque;
		}
	}

	double V = this->_config_info.box->V();
	_P = _T * (_N_constructs / V) + virial / (3. * V);
}

template<typename number>
string ConstructwisePressure<number>::get_output_string(llint curr_step) {
	update_pressure();
	return Utils::sformat("% .8e", get_P());
}

template class ConstructwisePressure<float> ;
template class ConstructwisePressure<double> ;

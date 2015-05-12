/*
 * Pressure.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#include "Pressure.h"

template<typename number>
Pressure<number>::Pressure(): BaseObservable<number>(), _stress_tensor(false) {

}

template<typename number>
Pressure<number>::~Pressure() {

}

template<typename number>
void Pressure<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	char raw_T[256];
	getInputString(&sim_inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number>(raw_T);

	int tmp = 0;
	getInputBoolAsInt(&my_inp, "stress_tensor", &tmp, 0);
	_stress_tensor = (bool) tmp;
}

template<typename number>
string Pressure<number>::get_output_string(llint curr_step) {
	number L = *this->_config_info.box_side;
	int N = *this->_config_info.N;

	std::vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();

	number virial = 0;
	number st_xx, st_yy, st_zz, st_xy, st_xz, st_yz;
	st_xx = st_yy = st_zz = st_xy = st_xz = st_yz = 0;
	// we loop on all the pairs in order to update the forces
	typename std::vector<ParticlePair<number> >::iterator it;
	for (it = pairs.begin(); it != pairs.end(); it ++ ) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;
		LR_vector<number> r = q->pos.minimum_image(p->pos, L);

		// pair_interaction will change these vectors, but we still need them in the next
		// first integration step. For this reason we copy and then restore their values
		// after the calculation
		LR_vector<number> old_p_force(p->force);
		LR_vector<number> old_q_force(q->force);
		LR_vector<number> old_p_torque(p->torque);
		LR_vector<number> old_q_torque(q->torque);

		p->force = LR_vector<number>(0, 0, 0);
		q->force = LR_vector<number>(0, 0, 0);
		p->torque = LR_vector<number>(0, 0, 0);
		q->torque = LR_vector<number>(0, 0, 0);

		this->_config_info.interaction->pair_interaction(p, q, NULL, true);

		st_xx += r.x * q->force.x;
		st_yy += r.y * q->force.y;
		st_zz += r.z * q->force.z;
		st_xy += r.x * q->force.y;
		st_xz += r.x * q->force.z;
		st_yz += r.y * q->force.z;

		virial += r*q->force;

		p->force = old_p_force;
		q->force = old_q_force;
		p->torque = old_p_torque;
		q->torque = old_q_torque;
	}

	number V = L*L*L;
	number P = N * _T / V + virial / (3*V);

	if(_stress_tensor) {
		for(int i = 0; i < *this->_config_info.N; i++) {
			BaseParticle<number> *p = this->_config_info.particles[i];
			st_xx += SQR(p->vel.x);
			st_yy += SQR(p->vel.y);
			st_zz += SQR(p->vel.z);
			st_xy += p->vel.x*p->vel.y;
			st_xz += p->vel.x*p->vel.z;
			st_yz += p->vel.y*p->vel.z;
		}

		return Utils::sformat("% 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf", P, st_xx / (3*V), st_yy / (3*V), st_zz / (3*V), st_xy / (3*V), st_xz / (3*V), st_yz / (3*V));
	}
	else return Utils::sformat("% 10.6lf", P);
}

template class Pressure<float>;
template class Pressure<double>;

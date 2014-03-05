/*
 * KineticEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "KineticEnergy.h"

template<typename number>
KineticEnergy<number>::KineticEnergy() {

}

template<typename number>
KineticEnergy<number>::~KineticEnergy() {

}

template<typename number>
number KineticEnergy<number>::get_kinetic_energy() {
	number K = (number) 0.f;
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(p->is_rigid_body()) K += (p->vel.norm() + p->L.norm()) * (number) 0.5f;
		else K += p->vel.norm() * (number) 0.5f;
	}
	K /= *this->_config_info.N;

	return K;
}

template<typename number>
std::string KineticEnergy<number>::get_output_string(llint curr_step) {
	number K = get_kinetic_energy();

	return Utils::sformat("% 10.6lf", K);
}

template class KineticEnergy<float>;
template class KineticEnergy<double>;

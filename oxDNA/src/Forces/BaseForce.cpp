/*
 * BaseForce.cpp
 *
 *  Created on: 22 giu 2017
 *      Author: lorenzo
 */

#include "BaseForce.h"
#include "../Particles/BaseParticle.h"

template <typename number>
BaseForce<number>::BaseForce () {
	_F0 = -1.;
	_rate = -1.;
	_direction = LR_vector<number>(1., 0., 0.);
	_pos0 = LR_vector<number>(0., 0., 0.);
	_site = -1;
	_stiff = 0.;
	_p_ptr = P_VIRTUAL;
}

template <typename number>
BaseForce<number>::~BaseForce () {

}

template <typename number>
void BaseForce<number>::_add_self_to_particles(BaseParticle<number> **particles, int N, std::string particle_string, std::string force_description) {
	std::vector<int> particle_indices_vector = Utils::getParticlesFromString(particles, N, particle_string, force_description.c_str());

	if (particle_indices_vector[0] != -1) {
		for (std::vector<int>::size_type i = 0; i < particle_indices_vector.size(); i++) {
			particles[particle_indices_vector[i]]->add_ext_force(this);
			OX_LOG (Logger::LOG_INFO, "Adding a %s on particle %d", force_description.c_str(), particle_indices_vector[i]);
		}
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding a %s on ALL particles", force_description.c_str());
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template class BaseForce<double>;
template class BaseForce<float>;

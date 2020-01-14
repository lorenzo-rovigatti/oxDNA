/*
 * BaseForce.cpp
 *
 *  Created on: 22 giu 2017
 *      Author: lorenzo
 */

#include "BaseForce.h"
#include "../Particles/BaseParticle.h"

BaseForce::BaseForce() {
	_F0 = -1.;
	_rate = -1.;
	_direction = LR_vector(1., 0., 0.);
	_pos0 = LR_vector(0., 0., 0.);
	_site = -1;
	_stiff = 0.;
	_p_ptr = P_VIRTUAL;
}

BaseForce::~BaseForce() {

}

void BaseForce::_add_self_to_particles(std::vector<BaseParticle *> &particles, int N, std::string particle_string, std::string force_description) {
	auto particle_indices_vector = Utils::getParticlesFromString(particles, N, particle_string, force_description.c_str());

	if(particle_indices_vector[0] != -1) {
		for(std::vector<int>::size_type i = 0; i < particle_indices_vector.size(); i++) {
			particles[particle_indices_vector[i]]->add_ext_force(ForcePtr(this));
			OX_LOG(Logger::LOG_INFO, "Adding a %s on particle %d", force_description.c_str(), particle_indices_vector[i]);
		}
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding a %s on ALL particles", force_description.c_str());
		for (int i = 0; i < N; i ++) {
			particles[i]->add_ext_force(ForcePtr(this));
		}
	}
}

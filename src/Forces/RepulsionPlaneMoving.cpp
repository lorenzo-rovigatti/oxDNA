/*
 * RepulsionPlaneMoving.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "RepulsionPlaneMoving.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"
#include "../Utilities/Utils.h"

#include <algorithm>

RepulsionPlaneMoving::RepulsionPlaneMoving() :
				BaseForce() {
	_particles_string = "-1";
	_ref_particles_string = "-1";
	_box_ptr = NULL;
	low_idx = high_idx = -1;
}

void RepulsionPlaneMoving::get_settings(input_file &inp) {
	getInputString(&inp, "particle", _particles_string, 1);
	getInputString(&inp, "ref_particle", _ref_particles_string, 1);

	getInputNumber(&inp, "stiff", &this->_stiff, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

void RepulsionPlaneMoving::init(BaseParticle ** particles, int N, BaseBox * box_ptr) {
	_box_ptr = box_ptr;

	auto particle_indices_vector = Utils::getParticlesFromString(particles, N, _particles_string, "moving repulsion plane force (RepulsionPlaneMoving.cpp)");
	auto ref_particle_indices_vector = Utils::getParticlesFromString(particles, N, _ref_particles_string, "moving repulsion plane force (RepulsionPlaneMoving.cpp)");

	sort(ref_particle_indices_vector.begin(), ref_particle_indices_vector.end());
	low_idx = ref_particle_indices_vector.front();
	high_idx = ref_particle_indices_vector.back();

	if((high_idx - low_idx + 1) != (int) ref_particle_indices_vector.size()) throw oxDNAException("RepulsionPlaneMoving requires the list of ref_particle indices to be contiguous");
	for(std::vector<int>::size_type i = 0; i < ref_particle_indices_vector.size(); i++) {
		ref_p_ptr.push_back(particles[ref_particle_indices_vector[i]]);
	}

	if(particle_indices_vector[0] != -1) {
		OX_LOG(Logger::LOG_INFO, "Adding repulsion_plane_moving force (RepulsionPlaneMoving.cpp) with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on particle %s", this->_stiff, this->_direction.x, this->_direction.y, this->_direction.z, _particles_string.c_str());
		for (std::vector<int>::size_type i = 0; i < particle_indices_vector.size(); i++) {
			particles[particle_indices_vector[i]]->add_ext_force(ForcePtr(this));
		}
	}
	else { // force affects all particles
		OX_LOG(Logger::LOG_INFO, "Adding repulsion_plane_moving force (RepulsionPlaneMoving.cpp) with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on ALL particles", this->_stiff, this->_direction.x, this->_direction.y, this->_direction.z);
		for (int i = 0; i < N; i ++) {
			particles[i]->add_ext_force(ForcePtr(this));
		}
	}
}

LR_vector RepulsionPlaneMoving::value(llint step, LR_vector &pos) {
	LR_vector force(0., 0., 0.);
	for(std::vector<int>::size_type i = 0; i < ref_p_ptr.size(); i++) {
		number distance = (pos - _box_ptr->get_abs_pos(ref_p_ptr[i])) * this->_direction;
		force += -this->_stiff * this->_direction * std::min(distance, (number) 0.);
	}
	return force;
}

number RepulsionPlaneMoving::potential(llint step, LR_vector &pos) {
	number V = 0.;
	for(std::vector<int>::size_type i = 0; i < ref_p_ptr.size(); i++) {
		number distance = (pos - _box_ptr->get_abs_pos(ref_p_ptr[i])) * this->_direction;
		V += (number) 0.5 * this->_stiff * distance * std::min(distance, (number) 0.);
	}
	return (number) V;
}

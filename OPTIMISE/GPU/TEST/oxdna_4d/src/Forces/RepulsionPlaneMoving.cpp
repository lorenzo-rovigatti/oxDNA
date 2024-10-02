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
	low_idx = high_idx = -1;
}

std::tuple<std::vector<int>, std::string> RepulsionPlaneMoving::init(input_file &inp) {
	BaseForce::init(inp);

	getInputString(&inp, "particle", _particles_string, 1);
	getInputString(&inp, "ref_particle", _ref_particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_direction.normalize();

	auto particle_indices_vector = Utils::get_particles_from_string(CONFIG_INFO->particles(), _particles_string, "moving repulsion plane force (RepulsionPlaneMoving.cpp)");
	auto ref_particle_indices_vector = Utils::get_particles_from_string(CONFIG_INFO->particles(), _ref_particles_string, "moving repulsion plane force (RepulsionPlaneMoving.cpp)");

	sort(ref_particle_indices_vector.begin(), ref_particle_indices_vector.end());
	low_idx = ref_particle_indices_vector.front();
	high_idx = ref_particle_indices_vector.back();

	if((high_idx - low_idx + 1) != (int) ref_particle_indices_vector.size()) throw oxDNAException("RepulsionPlaneMoving requires the list of ref_particle indices to be contiguous");
	for(std::vector<int>::size_type i = 0; i < ref_particle_indices_vector.size(); i++) {
		ref_p_ptr.push_back(CONFIG_INFO->particles()[ref_particle_indices_vector[i]]);
	}

	std::string description = Utils::sformat("RepulsionPlaneMoving force (RepulsionPlaneMoving.cpp) with stiffness %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]", _stiff, _direction.x, _direction.y, _direction.z);

	return std::make_tuple(particle_indices_vector, description);
}

LR_vector RepulsionPlaneMoving::value(llint step, LR_vector &pos) {
	LR_vector force(0., 0., 0.);
	for(std::vector<int>::size_type i = 0; i < ref_p_ptr.size(); i++) {
		number distance = (pos - CONFIG_INFO->box->get_abs_pos(ref_p_ptr[i])) * _direction;
		force += -_stiff * _direction * std::min(distance, (number) 0.);
	}
	return force;
}

number RepulsionPlaneMoving::potential(llint step, LR_vector &pos) {
	number V = 0.;
	for(std::vector<int>::size_type i = 0; i < ref_p_ptr.size(); i++) {
		number distance = (pos - CONFIG_INFO->box->get_abs_pos(ref_p_ptr[i])) * _direction;
		V += (number) 0.5 * _stiff * distance * std::min(distance, (number) 0.);
	}
	return (number) V;
}

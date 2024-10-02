/*
 * JqordanInteraction.cpp
 *
 *  Created on: 19/feb/2016
 *      Author: flavio
 */

#include "JordanInteraction.h"
#include "../Particles/JordanParticle.h"
#include "../Utilities/Utils.h"

JordanInteraction::JordanInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(JORDAN, _jordan_interaction);

	_s = (number) 0.3; // patch width
	_m = 6;            // LJ exponent
	_phi = M_PI / 6.;  // angle below the equator
	_int_k = 0.f;      // stiffness of the internal spring
	_SB_s = 1.;
	_n_patches = -1;

	_my_N = -1;
	_my_N3 = -1;
}

JordanInteraction::~JordanInteraction() {

}

void JordanInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	float tmp = 2.5;
	getInputFloat(&inp, "JORDAN_rcut", &tmp, 0);
	getInputInt(&inp, "JORDAN_N_patches", &_n_patches, 1);
	_rcut = (number) tmp;
	_sqr_rcut = _rcut * _rcut;

	getInputNumber(&inp, "JORDAN_phi", &_phi, 0);
	getInputNumber(&inp, "JORDAN_s", &_s, 0);
	getInputNumber(&inp, "JORDAN_int_k", &_int_k, 0);
	getInputNumber(&inp, "JORDAN_SB_s", &_SB_s, 0);
	getInputInt(&inp, "JORDAN_m", &_m, 0);
	//getInputString(&inp, "JORDAN_load_patches_from", patches_file_in, 0);
	//getInputString(&inp, "JORDAN_save_patches_to", patches_file_out, 1);
}

void JordanInteraction::init() {
	OX_LOG(Logger::LOG_INFO,"(JordanInteraction.cpp) Running Jordan interaction with s=%g, rcut=%g, m=%d, phi=%g, int_k=%g", _s, _rcut, _m, _phi, _int_k, _my_N3);
}

void JordanInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	OX_LOG(Logger::LOG_INFO,"(JordanInteraction.cpp) Allocating %d particles with 3 patches and %d with 4 patches", _my_N3, N - _my_N3);
	for(int i = 0; i < N; i++) {
		if(i > _my_N3 && _my_N3 >= 0)
			particles[i] = new JordanParticle(4, _phi, _int_k);
		else
			particles[i] = new JordanParticle(3, _phi, _int_k);
	}
}

number JordanInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number JordanInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number JordanInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _jordan_interaction(p, q, false, update_forces);
}

number JordanInteraction::_jordan_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {

	if(update_forces)
		throw oxDNAException("JordanInteraction cannot compute forces. Use MC/MC2/VMMC to simulate.");

	number rnorm = _computed_r.norm();

	// if outside of rcut, we return here
	if(rnorm > _sqr_rcut)
		return (number) 0.f;

	number rm = sqrt(rnorm);

	// otherwise, we keep going;
	number r_to_m = pow(rm, -(number) _m); // this returns (1./|r|)^m

	number energy = 4.f * (r_to_m * r_to_m - r_to_m);

	if(rm <= 1.f)
		return energy; // we do not consider the patches if r < 1

	// here we look for the closest patch on particle p;
	int imax1 = -1, imax2 = -1;
	number max1 = -2.f;
	for(uint i = 0; i < p->N_int_centers(); i++) {
		number tmp = ((p->int_centers[i] * _computed_r) / (0.5 * rm));
		if(tmp > max1) {
			max1 = tmp;
			imax1 = i;
		}
	}
	number max2 = -2.f;
	for(uint i = 0; i < q->N_int_centers(); i++) {
		number tmp = (-(q->int_centers[i] * _computed_r) / (0.5 * rm));
		if(tmp > max2) {
			max2 = tmp;
			imax2 = i;
		}
	}

	// now max1 and max2 are ready to be acos'd
	if(fabs(max1) >= 1.) {
		max1 = copysign(1., max1);
	}
	if(fabs(max2) >= 1.) {
		max2 = copysign(1., max2);
	}

	number theta1 = acos(max1);
	number theta2 = acos(max2);

	number fact = exp((-theta1 * theta1 - theta2 * theta2) / (2.f * _s * _s));

	// simmetry breaking term
	// we get a random vector and remove its projection onto r
	LR_vector pprime = p->int_centers[imax1] - _computed_r * ((p->int_centers[imax1] * _computed_r) / (rm * rm));
	LR_vector qprime = q->int_centers[imax2] - _computed_r * ((q->int_centers[imax2] * _computed_r) / (rm * rm));

	number Gamma_arg = (pprime * qprime) * (1.f / (pprime.module() * qprime.module()));
	if(fabs(Gamma_arg) > 1.f) {
		Gamma_arg = copysign(1.f, Gamma_arg);
	}
	number Gamma = acos(Gamma_arg);

	number fact_SB = exp(-Gamma * Gamma / (2.f * _SB_s * _SB_s));  // could be optimized joining with previous exponential

	return energy * fact * fact_SB;

}

void JordanInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);
	if(!topology.good()) {
		throw oxDNAException("(JordanInteraction.cpp) Can't read topology file '%s'. Aborting", _topology_filename);
	}
	std::string line;
	std::getline(topology, line);
	int check = sscanf(line.c_str(), "%d %d\n", &_my_N, &_my_N3);
	if(check != 2) {
		throw oxDNAException("(JordanInteraction.cpp) Wrong content in topology file '%s'. Expecting <N> and <N3>", _topology_filename);
	}
	topology.close();

	int N = particles.size();
	*N_strands = N;

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = i;
		particles[i]->type = 0;
		particles[i]->btype = 1;
		particles[i]->strand_id = i;
	}
}

void JordanInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

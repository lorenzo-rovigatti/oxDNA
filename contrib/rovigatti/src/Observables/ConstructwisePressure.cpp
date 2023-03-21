/*
 * ConstructwisePressure.cpp
 *
 *  Created on: 08/jan/2018
 *      Author: lorenzo
 */

#include "ConstructwisePressure.h"

#include "../Interactions/PolymerSwapInteraction.h"

ConstructwisePressure::ConstructwisePressure() :
				BaseObservable(),
				_P(0.),
				_shear_rate(0.),
				_construct_size(-1),
				_N_constructs(-1) {

}

ConstructwisePressure::~ConstructwisePressure() {

}

void ConstructwisePressure::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "construct_size", &_construct_size, 1);

	bool lees_edwards = false;
	getInputBool(&sim_inp, "lees_edwards", &lees_edwards, 0);
	if(lees_edwards) {
		getInputNumber(&sim_inp, "lees_edwards_shear_rate", &_shear_rate, 0);
	}

	std::string interaction;
	getInputString(&sim_inp, "interaction_type", interaction, 1);
	_is_polymer_swap = (interaction == "PolymerSwapInteraction");
}

void ConstructwisePressure::init() {
	BaseObservable::init();

	int N = _config_info->N();

	if(N % _construct_size) {
		throw oxDNAException("ConstructwisePressure: the total number of particles (%d) is not a multiple of the construct size specified in the input file (%d)", N, _construct_size);
	}
	_N_constructs = N / _construct_size;
	_construct_coms.resize(_N_constructs);

	if(_is_polymer_swap) {
		OX_LOG(Logger::LOG_INFO, "ConstructwisePressure: PolymerSwapInteraction detected");
	}
}

void ConstructwisePressure::update_pressure() {
	fill(_construct_coms.begin(), _construct_coms.end(), LR_vector());
	for(auto p: _config_info->particles()) {
		int p_construct = p->index / _construct_size;
		_construct_coms[p_construct] += _config_info->box->get_abs_pos(p);
	}
	for(int i = 0; i < _N_constructs; i++) {
		_construct_coms[i] /= _construct_size;
	}

	std::vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();

	double virial = 0;
	double energy = 0.;
	// we loop over all the pairs in order to update the forces
	for(auto &pair: pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;

		int p_construct = p->index / _construct_size;
		int q_construct = q->index / _construct_size;

		if(p_construct != q_construct) {
			LR_vector r = _config_info->box->min_image(_construct_coms[p_construct], _construct_coms[q_construct]);

			// pair_interaction will change these vectors, but we still need them in the next
			// first integration step. For this reason we copy and then restore their values
			// after the calculation
			LR_vector old_p_force(p->force);
			LR_vector old_q_force(q->force);
			LR_vector old_p_torque(p->torque);
			LR_vector old_q_torque(q->torque);

			p->force = q->force = p->torque = q->torque = LR_vector();

			energy += (double) _config_info->interaction->pair_interaction(p, q, true, true);

			virial -= (r * p->force);

			p->force = old_p_force;
			q->force = old_q_force;
			p->torque = old_p_torque;
			q->torque = old_q_torque;
		}
	}

	double V = _config_info->box->V();
	_P = _config_info->temperature() * (_N_constructs / V) + virial / (3. * V);
}

void ConstructwisePressure::update_pressure_PolymerSwap() {
	PolymerSwapInteraction *interaction = dynamic_cast<PolymerSwapInteraction *>(_config_info->interaction);
	interaction->begin_energy_and_force_computation();

	std::vector<LR_vector> forces;
	forces.reserve(CONFIG_INFO->N());
	auto particles = CONFIG_INFO->particles();
	// backup the particles' forces
	std::transform(particles.begin(), particles.end(), std::back_inserter(forces), [](BaseParticle *p) { return p->force; });
	
	std::vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();

	// we loop over all the pairs in order to update the virial stored in the interaction
	for(auto &pair : pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;
		interaction->pair_interaction(p, q, true, true);
	}

	// copy back the backupped forces
	for(auto p : particles) {
		p->force = forces[p->index];
	}

	_P = interaction->P_inter_chain();
}

std::string ConstructwisePressure::get_output_string(llint curr_step) {
	if(_is_polymer_swap) {
		update_pressure_PolymerSwap();
	}
	else {
		update_pressure();
	}
	return Utils::sformat("% .8e", get_P());
}


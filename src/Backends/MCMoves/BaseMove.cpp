/*
 * BaseMove.cpp
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#include "BaseMove.h"

BaseMove::BaseMove() :
				_Info(ConfigInfo::instance()) {
	_attempted = 0;
	_accepted = 0;
	_T = -1.;
	prob = (number) 1.f;
	_target_acc_rate = 0.25;
	_equilibration_steps = 0;
	_adjust_moves = false;
	_compute_energy_before = true;
	_restrict_to_type = -1;
	_rej_fact = 1.001;
	_acc_fact = 0.;

	CONFIG_INFO->subscribe("T_updated", [this]() { this->_on_T_update(); });
}

BaseMove::~BaseMove() {

}

void BaseMove::_on_T_update() {
	_T = _Info->temperature();
}

void BaseMove::get_settings(input_file &inp, input_file &sim_inp) {
	getInputString(&inp, "type", _name, 1);
	getInputNumber(&inp, "prob", &prob, 0);
	getInputLLInt(&sim_inp, "equilibration_steps", &_equilibration_steps, 0);
	getInputBool(&inp, "adjust_moves", &_adjust_moves, 0);
	getInputNumber(&inp, "target_acc_rate", &_target_acc_rate, 0);
	getInputBool(&inp, "compute_energy_before", &_compute_energy_before, 0);
	getInputInt(&inp, "restrict_to_type", &_restrict_to_type, 0);
}

void BaseMove::init() {
	_T = _Info->temperature();
	_acc_fact = 1. + 0.001 * (_target_acc_rate - 1.) / (0.001 - _target_acc_rate);
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h) BaseMove init for move type %s;", _name.c_str());
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h)               (b, a) = (%g, %g), target_acc_rate=%g, adjust_moves=%d,", (_rej_fact - 1.), (_acc_fact - 1.), _target_acc_rate, _adjust_moves);
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h)               compute_energy_before=%d, restrict_to_type=%d", _compute_energy_before, _restrict_to_type);
}

void BaseMove::log_parameters() {
	OX_LOG(Logger::LOG_INFO, "Logging parameters for MC2 Move with type = %s", _name.c_str());
	OX_LOG(Logger::LOG_INFO, "\tprob = %g", prob);
	OX_LOG(Logger::LOG_INFO, "\tadjust_moves = %d", int(_adjust_moves));
	OX_LOG(Logger::LOG_INFO, "\ttarget_acc_rate = %g", _target_acc_rate);
	OX_LOG(Logger::LOG_INFO, "\tcompute_energy_before = %d", int(_compute_energy_before));
	OX_LOG(Logger::LOG_INFO, "\trestrict_to_type = %d", int(_restrict_to_type));
}

number BaseMove::particle_energy(BaseParticle * p) {
	number res = (number) 0.f;
	typename vector<ParticlePair>::iterator it = p->affected.begin();
	for(; it != p->affected.end(); it++) {
		number de = _Info->interaction->pair_interaction_bonded(it->first, it->second);
		res += de;
	}
	if(_Info->interaction->get_is_infinite() == true){
		return (number) 1.e12;
	}

	std::vector<BaseParticle *> neighs = _Info->lists->get_neigh_list(p);
	for(unsigned int n = 0; n < neighs.size(); n++) {
		BaseParticle *q = neighs[n];
		res += _Info->interaction->pair_interaction_nonbonded(p, q);
		if(_Info->interaction->get_is_infinite() == true) {
			return (number) 1.e12;
		}
	}
	return res;
}

number BaseMove::system_energy() {
	number res = (number) 0.f;

	for(auto p: _Info->particles()) {
		if(p->n3 != P_VIRTUAL) res += _Info->interaction->pair_interaction_bonded(p, p->n3);
		// we omit E(p,p->n5) because it gets counted as E(q, q->n3);
		std::vector<BaseParticle *> neighs = _Info->lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			if(p->index < q->index) {
				res += _Info->interaction->pair_interaction_nonbonded(p, q);
				if(_Info->interaction->get_is_infinite() == true) {
					return (number) 1.e12;
				}
			}
		}
	}

	return res;
}

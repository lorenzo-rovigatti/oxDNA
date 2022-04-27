/*
 * StretchedBonds.cpp
 *
 *  Created on: Mar 20, 2016
 *      Author: Flavio
 */

#include "StretchedBonds.h"
#include "../Interactions/DNAInteraction.h"

StretchedBonds::StretchedBonds() {
	_print_list = false;
	_threshold = 1.f;
}

StretchedBonds::~StretchedBonds() {

}

void StretchedBonds::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputBool(&my_inp, "print_list", &_print_list, 0);
	getInputNumber(&my_inp, "threshold", &_threshold, 0);
}

void StretchedBonds::init() {
	BaseObservable::init();
	OX_LOG(Logger::LOG_INFO, "(StretchedBonds.cpp) Initialized StretchedBonds observables with print_list = %d and threshold = %g", (int)_print_list, _threshold);
	return;
}

std::string StretchedBonds::get_output_string(llint curr_step) {
	std::string ret;
	std::vector<int> stretched;
	int num = 0;
	int N = _config_info->N();
	for(int i = 0; i < N; i++) {
		BaseParticle * p = _config_info->particles()[i];
		number energy;
		if(p->n3 == P_VIRTUAL) continue;

		energy = _config_info->interaction->pair_interaction_term(DNAInteraction::BACKBONE, p, p->n3);

		if(energy > _threshold) {
			num++;
			stretched.push_back(p->index);
			stretched.push_back(p->n3->index);
		}
	}

	ret = Utils::sformat("%d", num);

	if(_print_list) {
		for(int i = 0; i < (int) stretched.size(); i += 2) {
			ret += Utils::sformat(" %d-%d", stretched[i], stretched[i + 1]);
		}
	}

	return ret;

}

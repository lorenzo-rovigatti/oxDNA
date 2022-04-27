/*
 * Gyradius.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include "Gyradius.h"
#include "Utilities/Utils.h"

Gyradius::Gyradius() {
	_accumulate = false;
	_avg_gr2 = 0;
	_counter = 0;
	_type = -1;
}

Gyradius::~Gyradius() {

}

void Gyradius::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	int tmp = 0;
	getInputBoolAsInt(&my_inp, "accumulate", &tmp, 0);
	_accumulate = (bool) tmp;

	getInputInt(&my_inp, "rg_type", &_type, 0);
	if(_type != -1 && _type != P_A && _type != P_B) throw oxDNAException("rg_type should be either -1 (all particles), 0 (P_A) or 1 (P_B)");
}

std::string Gyradius::get_output_string(llint curr_step) {
	int N = _config_info->N();

	number rg2 = 0;
	int N_type = 0;
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		if(p->type != 2 && (_type == -1 || p->type == _type)) {
			N_type++;
			for(int j = i + 1; j < N; j++) {
				BaseParticle *q = _config_info->particles()[j];
				if(_type == -1 || q->type == _type) rg2 += _config_info->box->sqr_min_image_distance(q->pos, p->pos);
			}
		}
	}
	rg2 /= SQR(N_type);

	if(_accumulate) {
		_avg_gr2 += rg2;
		_counter++;
		return Utils::sformat("%f %d", sqrt(_avg_gr2 / _counter), N_type);
	}

	return Utils::sformat("%f %d", sqrt(rg2), N_type);
}

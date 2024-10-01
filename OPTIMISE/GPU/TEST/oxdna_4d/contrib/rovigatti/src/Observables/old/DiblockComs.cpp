/*
 * DiblockComs.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include <sstream>

#include "DiblockComs.h"
#include "Utilities/Utils.h"

DiblockComs::DiblockComs() {
	_only_intra = false;
}

DiblockComs::~DiblockComs() {

}

void DiblockComs::get_settings(input_file &my_inp, input_file &sim_inp) {
	int only_intra = 0;
	getInputBoolAsInt(&my_inp, "only_intra", &only_intra, 0);
	_only_intra = (bool) only_intra;

	CHECK_BOX("DiblockComs", my_inp);
}

void DiblockComs::init() {
	BaseObservable::init();
}

std::string DiblockComs::get_output_string(llint curr_step) {
	int N = _config_info->N();

	LR_vector coms[2][2];
	int counters[2][2] = { { 0, 0 }, { 0, 0 } };

	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		if(p->strand_id > 1) throw oxDNAException("The system should contain just two chains");
		if(p->type != P_A && p->type != P_B) throw oxDNAException("Only A and B particle types allowed");
		coms[p->strand_id][p->type] += _config_info->box->get_abs_pos(p);
		counters[p->strand_id][p->type]++;
	}
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			coms[i][j] /= counters[i][j];

	std::stringstream ret;
	if(!_only_intra) {
		number r_AA = sqrt(_config_info->box->sqr_min_image_distance(coms[1][P_A], coms[0][P_A]));
		number r_AB = sqrt(_config_info->box->sqr_min_image_distance(coms[1][P_B], coms[0][P_A]));
		number r_BA = sqrt(_config_info->box->sqr_min_image_distance(coms[1][P_A], coms[0][P_B]));
		number r_BB = sqrt(_config_info->box->sqr_min_image_distance(coms[1][P_B], coms[0][P_B]));

		ret << r_AA << " " << r_AB << " " << r_BA << " " << r_BB << " ";
	}

	number r_intra = sqrt(_config_info->box->sqr_min_image_distance(coms[0][P_B], coms[0][P_A]));
	ret << r_intra;

	return ret.str();
}

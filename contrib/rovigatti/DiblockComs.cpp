/*
 * DiblockComs.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include <sstream>

#include "DiblockComs.h"
#include "../Utilities/Utils.h"

template<typename number>
DiblockComs<number>::DiblockComs() {
	_only_intra = false;
}

template<typename number>
DiblockComs<number>::~DiblockComs() {

}

template<typename number>
void DiblockComs<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int only_intra = 0;
	getInputBoolAsInt(&my_inp, "only_intra", &only_intra, 0);
	_only_intra = (bool) only_intra;
}

template<typename number>
void DiblockComs<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);
}

template<typename number>
std::string DiblockComs<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	number L = *this->_config_info.box_side;

	LR_vector<number> coms[2][2];
	int counters[2][2] = {{0, 0}, {0, 0}};

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(p->strand_id > 1) throw oxDNAException("The system should contain just two chains");
		if(p->type != P_A && p->type != P_B) throw oxDNAException("Only A and B particle types allowed");
		coms[p->strand_id][p->type] += p->get_abs_pos(L);
		counters[p->strand_id][p->type]++;
	}
	for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) coms[i][j] /= counters[i][j];

	stringstream ret;
	if(!_only_intra) {
		number r_AA = sqrt(coms[0][P_A].sqr_min_image_distance(coms[1][P_A], L));
		number r_AB = sqrt(coms[0][P_A].sqr_min_image_distance(coms[1][P_B], L));
		number r_BA = sqrt(coms[0][P_B].sqr_min_image_distance(coms[1][P_A], L));
		number r_BB = sqrt(coms[0][P_B].sqr_min_image_distance(coms[1][P_B], L));

		ret << r_AA << " " << r_AB << " " << r_BA << " " << r_BB << " ";
	}

	number r_intra = sqrt(coms[0][P_A].sqr_min_image_distance(coms[0][P_B], L));
	ret << r_intra;

	return ret.str();
}

template class DiblockComs<float>;
template class DiblockComs<double>;

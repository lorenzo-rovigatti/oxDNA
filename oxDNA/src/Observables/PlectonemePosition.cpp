/*
 * DenaturationPattern.cpp
 *
 *  Created on: Aug 19, 2013
 *      Author: matek
 */

#include "PlectonemePosition.h"
#include "../Interactions/DNAInteraction.h"
#include "../Utilities/OrderParameters.h"

template<typename number>
PlectonemePosition<number>::PlectonemePosition() {
	_critical_bp_distance = 40;
	_critical_strand_distance = 8.5;
	_midpoints = NULL;
}

template<typename number>
PlectonemePosition<number>::~PlectonemePosition() {
	if(_midpoints != NULL) delete[] _midpoints;
}

template<typename number>
void PlectonemePosition<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputInt(&my_inp, "bp_dist_along_str", &_critical_bp_distance, 0);
	getInputNumber(&my_inp, "dist_bet_str", &_critical_strand_distance, 0);
}

template<typename number> void
PlectonemePosition<number>::init(ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);

   _midpoints = new LR_vector<number>[*this->_config_info.N / 2];
}

template<typename number>
std::string PlectonemePosition<number>::get_output_string(llint curr_step) {

	int max_size = -1, max_plectoneme_start = -1, max_plectoneme_end = -1;

	for(int i = 0; i < *this->_config_info.N / 2; i++) {
		_midpoints[i] = (this->_config_info.particles[i]->pos + this->_config_info.particles[*this->_config_info.N - 1 - i]->pos);
	}

	for(int i = 0; i < *this->_config_info.N / 2 - _critical_bp_distance; ++i) {
		int cur_plecto_start = -1, cur_plecto_end = -1;
		for(int j = i + _critical_bp_distance; j < *this->_config_info.N / 2; ++j) {
			LR_vector<number> cur_distvec = _midpoints[j] - _midpoints[i];
			if(cur_distvec.module() < _critical_strand_distance) {

				if(cur_plecto_start == -1)
					cur_plecto_start = i;
			}
			else {
				if(cur_plecto_start != -1) {
					cur_plecto_end = j - 1;
					if(max_size < (cur_plecto_end - cur_plecto_start)) {
						max_size = cur_plecto_end - cur_plecto_start;
						max_plectoneme_start = cur_plecto_start;
						max_plectoneme_end = cur_plecto_end;
					}
					i = j - 1;
					break;
				}
			}
		}
	}

	//Format output string
	std::string plectopos_string;
	char outstr[100];
	if(max_plectoneme_start != -1 && max_plectoneme_end != -1) {
		number plecto_pos = 0.5 * ((number) max_plectoneme_start + (number) max_plectoneme_end);
		number plecto_size = (number) max_plectoneme_end - (number) max_plectoneme_start;
		sprintf(outstr, "%f %d %f %f", _critical_strand_distance, _critical_bp_distance, plecto_pos, plecto_size);
		plectopos_string = outstr;
	}
	else sprintf(outstr, "%f %d nan nan", _critical_strand_distance, _critical_bp_distance);
	plectopos_string = outstr;

	return plectopos_string;

}

template class PlectonemePosition<float> ;
template class PlectonemePosition<double> ;

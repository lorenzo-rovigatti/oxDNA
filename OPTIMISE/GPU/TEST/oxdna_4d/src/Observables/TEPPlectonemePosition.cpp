/*
 * TEPPlectonemePosition.cpp
 *
 *  Created on: June 21, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "TEPPlectonemePosition.h"

#include "../Interactions/DNAInteraction.h"
#include "../Utilities/OrderParameters.h"

#include <sstream>

TEPPlectonemePosition::TEPPlectonemePosition() {
	_bead_minimum_distance = 7;
	_distance_threshold = 2.;
	_old_pos = -1.;
	_print_pos = false;
}

TEPPlectonemePosition::~TEPPlectonemePosition() {

}

void TEPPlectonemePosition::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "bead_minimum_distance", &_bead_minimum_distance, 0);
	getInputNumber(&my_inp, "distance_threshold", &_distance_threshold, 0);
	getInputBool(&my_inp, "print_position", &_print_pos, 0);
}

std::string TEPPlectonemePosition::get_output_string(llint curr_step) {
	// going from left to right
	int left_max_size = -1;
	int left_max_plectoneme_start = -1;
	int left_max_plectoneme_end = -1;
	for(int i = 0; i < _config_info->N() - _bead_minimum_distance; ++i) {
		int cur_plecto_start = -1;
		int cur_plecto_end = -1;
		for(int j = i + _bead_minimum_distance; j < _config_info->N(); ++j) {
			LR_vector cur_distvec = _config_info->particles()[j]->pos - _config_info->particles()[i]->pos;
			if(cur_distvec.module() < _distance_threshold) {
				if(cur_plecto_start == -1) cur_plecto_start = i;
			}
			else {
				if(cur_plecto_start != -1) {
					cur_plecto_end = j - 1;
					if(left_max_size < (cur_plecto_end - cur_plecto_start)) {
						left_max_size = cur_plecto_end - cur_plecto_start;
						left_max_plectoneme_start = cur_plecto_start;
						left_max_plectoneme_end = cur_plecto_end;
					}
					i = j - 1;
					break;
				}
			}
		}
	}

	// going from right to left
	int right_max_size = -1;
	int right_max_plectoneme_start = -1;
	int right_max_plectoneme_end = -1;
	for(int i = _config_info->N() - _bead_minimum_distance; i >= 0; --i) {
		int cur_plecto_start = -1;
		int cur_plecto_end = -1;
		for(int j = i - _bead_minimum_distance; j >= 0; --j) {
			LR_vector cur_distvec = _config_info->particles()[j]->pos - _config_info->particles()[i]->pos;
			if(cur_distvec.module() < _distance_threshold) {
				if(cur_plecto_start == -1) cur_plecto_start = i;
			}
			else {
				if(cur_plecto_start != -1) {
					cur_plecto_end = j + 1;
					if(right_max_size < (cur_plecto_start - cur_plecto_end)) {
						right_max_size = cur_plecto_start - cur_plecto_end;
						// we define start and end so that end - start > 0
						right_max_plectoneme_start = cur_plecto_end;
						right_max_plectoneme_end = cur_plecto_start;
					}
					i = j + 1;
					break;
				}
			}
		}
	}

	stringstream ss;
	number plecto_pos = -1;
	number left_plecto_pos = -1;
	number right_plecto_pos = -1;
	number plecto_size = -1;
	if(left_max_plectoneme_start != -1 && left_max_plectoneme_end != -1 && right_max_plectoneme_start != -1 && right_max_plectoneme_end != -1) {
		left_plecto_pos = 0.5 * ((number) left_max_plectoneme_start + (number) left_max_plectoneme_end);
		number left_plecto_size = (number) left_max_plectoneme_end - (number) left_max_plectoneme_start;
		right_plecto_pos = 0.5 * ((number) right_max_plectoneme_start + (number) right_max_plectoneme_end);
		number right_plecto_size = (number) right_max_plectoneme_end - (number) right_max_plectoneme_start;
		if(_old_pos < 0.) {
			plecto_pos = (left_plecto_pos + right_plecto_pos) / 2.;
			plecto_size = (left_plecto_size + right_plecto_size) / 2.;
		}
		else {
			int delta_left = abs(left_plecto_pos - _old_pos);
			int delta_right = abs(right_plecto_pos - _old_pos);
			if(delta_left < delta_right) {
				plecto_pos = left_plecto_pos;
				plecto_size = left_plecto_size;
			}
			else {
				plecto_pos = right_plecto_pos;
				plecto_size = right_plecto_size;
			}
		}
		_old_pos = plecto_pos;
	}
	ss.precision(6);
	ss << plecto_pos << " " << left_plecto_pos << " " << right_plecto_pos << " " << plecto_size;
	if(_print_pos) {
		int idx = round(plecto_pos);
		LR_vector pos(0., 0., 0.);
		if(idx > 0 && idx < _config_info->N()) {
			BaseParticle *p = _config_info->particles()[idx];
			pos = p->pos;
		}
		ss << " " << pos.x << " " << pos.y << " " << pos.z;
	}
	return ss.str();

}

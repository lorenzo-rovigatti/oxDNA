/*
 * Pitch.cpp
 *
 *  Created on: Mar 14, 2014
 *      Author: Ben Snodin
 */

#include "Pitch.h"
#include "../Particles/DNANucleotide.h"

#include <sstream>

template<typename number>
Pitch<number>::Pitch() {
}

template<typename number>
Pitch<number>::~Pitch() {
}

template<typename number>
void Pitch<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int tmp = 0;

	tmp = 0;
	getInputInt(&my_inp,"bp1a_id", &tmp, 1);
	_bp1a_id = tmp;

	tmp = 0;
	getInputInt(&my_inp,"bp1b_id", &tmp, 1);
	_bp1b_id = tmp;

	tmp = 0;
	getInputInt(&my_inp,"bp2a_id", &tmp, 1);
	_bp2a_id = tmp;

	tmp = 0;
	getInputInt(&my_inp,"bp2b_id", &tmp, 1);
	_bp2b_id = tmp;
}

template<typename number>
std::string Pitch<number>::get_output_string(llint curr_step) {
	std::stringstream output_str;

	BaseParticle<number> *bp1a = this->_config_info.particles[_bp1a_id];
	BaseParticle<number> *bp1b = this->_config_info.particles[_bp1b_id];
	BaseParticle<number> *bp2a = this->_config_info.particles[_bp2a_id];
	BaseParticle<number> *bp2b = this->_config_info.particles[_bp2b_id];

	// base-base vector for each base pair; vector from a to b
	LR_vector<number> bp1_rbase = bp1b->pos.minimum_image(bp1a->pos, *(this->_config_info.box_side)) + bp1b->int_centers[DNANucleotide<number>::BASE] - bp1a->int_centers[DNANucleotide<number>::BASE];
	LR_vector<number> bp2_rbase = bp2b->pos.minimum_image(bp2a->pos, *(this->_config_info.box_side)) + bp2b->int_centers[DNANucleotide<number>::BASE] - bp2a->int_centers[DNANucleotide<number>::BASE];

	// base-base midpoint for each base pair; vector from 1a to 2a and 1b to 2b
	LR_vector<number> bpa_r = bp2a->pos.minimum_image(bp1a->pos, *(this->_config_info.box_side)) + bp2a->int_centers[DNANucleotide<number>::BASE] - bp1a->int_centers[DNANucleotide<number>::BASE];
	LR_vector<number> bpb_r = bp2b->pos.minimum_image(bp1b->pos, *(this->_config_info.box_side)) + bp2b->int_centers[DNANucleotide<number>::BASE] - bp1b->int_centers[DNANucleotide<number>::BASE];
	
	LR_vector<number> helix_axis = (bpa_r + bpb_r) * 0.5;
	number helix_axis_mod = helix_axis.module();
	LR_vector<number> helix_axis_dir = helix_axis / helix_axis_mod;
	
	// project base-base vectors into plane perpendicular to helix axis
	LR_vector<number> bp1_rbase_projected = bp1_rbase - helix_axis_dir * (bp1_rbase * helix_axis_dir);
	LR_vector<number> bp2_rbase_projected = bp2_rbase - helix_axis_dir * (bp2_rbase * helix_axis_dir);

	LR_vector<number> bp1_rbase_projected_dir = bp1_rbase_projected / bp1_rbase_projected.module();
	LR_vector<number> bp2_rbase_projected_dir = bp2_rbase_projected / bp2_rbase_projected.module();

	// get angle between projected base-base vectors
	number dot_product = bp1_rbase_projected_dir * bp2_rbase_projected_dir;
	number theta = acos(dot_product);

	output_str << theta;
	return output_str.str();
}

template class Pitch<float>;
template class Pitch<double>;

/*
 * Pitch.cpp
 *
 *  Created on: Mar 14, 2014
 *      Author: Ben Snodin
 */

#include "Pitch.h"
#include "../Particles/DNANucleotide.h"

#include <sstream>

Pitch::Pitch() {
}

Pitch::~Pitch() {
}

void Pitch::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	int tmp = 0;

	tmp = 0;
	getInputInt(&my_inp, "bp1a_id", &tmp, 1);
	_bp1a_id = tmp;

	tmp = 0;
	getInputInt(&my_inp, "bp1b_id", &tmp, 1);
	_bp1b_id = tmp;

	tmp = 0;
	getInputInt(&my_inp, "bp2a_id", &tmp, 1);
	_bp2a_id = tmp;

	tmp = 0;
	getInputInt(&my_inp, "bp2b_id", &tmp, 1);
	_bp2b_id = tmp;
}

std::string Pitch::get_output_string(llint curr_step) {
	std::stringstream output_str;

	BaseParticle *bp1a = _config_info->particles()[_bp1a_id];
	BaseParticle *bp1b = _config_info->particles()[_bp1b_id];
	BaseParticle *bp2a = _config_info->particles()[_bp2a_id];
	BaseParticle *bp2b = _config_info->particles()[_bp2b_id];

	BaseBox * mybox = _config_info->box;

	// base-base vector for each base pair; vector from a to b
	LR_vector bp1_rbase = mybox->min_image(bp1a->pos, bp1b->pos) + bp1b->int_centers[DNANucleotide::BASE] - bp1a->int_centers[DNANucleotide::BASE];
	LR_vector bp2_rbase = mybox->min_image(bp2a->pos, bp2b->pos) + bp2b->int_centers[DNANucleotide::BASE] - bp2a->int_centers[DNANucleotide::BASE];

	// base-base midpoint for each base pair; vector from 1a to 2a and 1b to 2b
	LR_vector bpa_r = mybox->min_image(bp1a->pos, bp2a->pos) + bp2a->int_centers[DNANucleotide::BASE] - bp1a->int_centers[DNANucleotide::BASE];
	LR_vector bpb_r = mybox->min_image(bp1b->pos, bp2b->pos) + bp2b->int_centers[DNANucleotide::BASE] - bp1b->int_centers[DNANucleotide::BASE];

	LR_vector helix_axis = (bpa_r + bpb_r) * 0.5;
	number helix_axis_mod = helix_axis.module();
	LR_vector helix_axis_dir = helix_axis / helix_axis_mod;

	// project base-base vectors into plane perpendicular to helix axis
	LR_vector bp1_rbase_projected = bp1_rbase - helix_axis_dir * (bp1_rbase * helix_axis_dir);
	LR_vector bp2_rbase_projected = bp2_rbase - helix_axis_dir * (bp2_rbase * helix_axis_dir);

	LR_vector bp1_rbase_projected_dir = bp1_rbase_projected / bp1_rbase_projected.module();
	LR_vector bp2_rbase_projected_dir = bp2_rbase_projected / bp2_rbase_projected.module();

	// get angle between projected base-base vectors
	number dot_product = bp1_rbase_projected_dir * bp2_rbase_projected_dir;
	number theta = acos(dot_product);

	output_str << theta;
	return output_str.str();
}

/*
 * meta_utils.h
 *
 *  Created on: Apr 3, 2022
 *      Author: lorenzo
 */

#ifndef SRC_FORCES_METADYNAMICS_META_UTILS_H_
#define SRC_FORCES_METADYNAMICS_META_UTILS_H_

#include "../../Particles/BaseParticle.h"
#include "../../Utilities/Utils.h"

#include <numeric>

namespace meta {

inline number interpolate_potential(const number my_x, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	int ix_left = static_cast<int>(std::floor((my_x - xmin) / dX));
	int ix_right = ix_left + 1;

    number x_left  = xmin + ix_left * dX;

    number w_right = (my_x - x_left) / dX;
    number w_left  = 1.0 - w_right;

    return w_left * potential_grid[ix_left] + w_right * potential_grid[ix_right];
}

inline number get_x_force(const number my_x, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	int ix_left = static_cast<int>(std::floor((my_x - xmin) / dX));
	int ix_right = ix_left + 1;
	return -(potential_grid[ix_right] - potential_grid[ix_left]) / dX;
}

LR_vector particle_list_com(const std::vector<BaseParticle *> &list, BaseBox *box_ptr);

std::tuple<std::vector<int>, std::vector<BaseParticle *>> get_particle_lists(input_file &inp, std::string key, std::vector<BaseParticle *> &particles, std::string description);

std::vector<number> split_to_numbers(const std::string &str, const std::string &delims);

struct CoordSettings {
	int coord_mode;
	number hb_energy_cutoff, hb_transition_width, d0, r0;
	int n;

	CoordSettings();

	void get_settings(input_file &inp);
};

LR_vector distance(std::pair<BaseParticle *, BaseParticle *> &pair);
number coordination(CoordSettings &settings, std::vector<std::pair<BaseParticle *, BaseParticle *>> &all_pairs);
number get_pair_contribution(CoordSettings &settings, std::pair<BaseParticle*, BaseParticle*> &pair);
number smooth_hb_contribution(number hb_energy_cutoff, number hb_transition_width, number hb_energy);

}

#endif /* SRC_FORCES_METADYNAMICS_META_UTILS_H_ */

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

}

#endif /* SRC_FORCES_METADYNAMICS_META_UTILS_H_ */

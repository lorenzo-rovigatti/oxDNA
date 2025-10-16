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
	number x_left = dX * std::floor(my_x / dX);
	number x_right = x_left + dX;
	number ix_left = std::floor((my_x - xmin) / dX);
	number ix_right = ix_left + 1;
	number f11 = potential_grid[ix_left];
	number f21 = potential_grid[ix_right];
	number fx = (x_right - my_x) / dX * f11 + (my_x - x_left) / dX * f21;
	return fx;
}

inline number get_x_force(const number x, const number dX, const number xmin, const std::vector<number> &potential_grid) {
	number ix_left = std::floor((x - xmin) / dX);
	number ix_right = ix_left + 1;
	return -(potential_grid[ix_right] - potential_grid[ix_left]) / dX;
}

LR_vector particle_list_com(const std::vector<BaseParticle *> &list, BaseBox *box_ptr);

std::tuple<std::vector<int>, std::vector<BaseParticle *>> get_particle_lists(input_file &inp, std::string key, std::vector<BaseParticle *> &particles, std::string description);

std::vector<number> split_to_numbers(const std::string &str, const std::string &delims);

}

#endif /* SRC_FORCES_METADYNAMICS_META_UTILS_H_ */

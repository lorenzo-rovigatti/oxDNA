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

LR_vector particle_list_com(const std::vector<BaseParticle *> &list, BaseBox *box_ptr);

std::tuple<std::vector<int>, std::vector<BaseParticle *>> get_particle_lists(input_file &inp, std::string key, std::vector<BaseParticle *> &particles, std::string description);

std::vector<number> split_to_numbers(const std::string &str, const std::string &delims);

}

#endif /* SRC_FORCES_METADYNAMICS_META_UTILS_H_ */

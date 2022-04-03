/*
 * meta_utils.h
 *
 *  Created on: Apr 3, 2022
 *      Author: lorenzo
 */

#ifndef SRC_FORCES_METADYNAMICS_META_UTILS_H_
#define SRC_FORCES_METADYNAMICS_META_UTILS_H_

#include "../../Particles/BaseParticle.h"

#include <numeric>

LR_vector particle_list_com(std::vector<BaseParticle *> list, BaseBox *box_ptr) {
	auto sum_abs_pos = [&](LR_vector sum, BaseParticle *p) { return sum + box_ptr->get_abs_pos(p); };
	return std::accumulate(list.begin(), list.end(), LR_vector(), sum_abs_pos) / (number) list.size();
}



#endif /* SRC_FORCES_METADYNAMICS_META_UTILS_H_ */

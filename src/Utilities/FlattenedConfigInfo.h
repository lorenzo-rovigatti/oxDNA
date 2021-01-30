/*
 * FlattenedConfigInfo.h
 *
 *  Created on: Jan 30, 2021
 *      Author: lorenzo
 */

#ifndef SRC_UTILITIES_FLATTENEDCONFIGINFO_H_
#define SRC_UTILITIES_FLATTENEDCONFIGINFO_H_

#include "../defs.h"

#include <vector>
#include <memory>

class BaseParticle;

struct FlattenedVectorArray {
	int rows();
	int cols();
	size_t size();

	void set(int p_idx, LR_vector &v);

	std::vector<number> data;
};

struct FlattenedConfigInfo {
	void update(long long int step, const std::vector<BaseParticle *> &particles);

	FlattenedVectorArray positions;
	FlattenedVectorArray a1s, a3s;
	std::vector<int> types;

	long long int last_updated = -1;
};

#endif /* SRC_UTILITIES_FLATTENEDCONFIGINFO_H_ */

/*
 * FlattenedConfigInfo.cpp
 *
 *  Created on: Jan 30, 2021
 *      Author: lorenzo
 */

#include "FlattenedConfigInfo.h"

#include "../Particles/BaseParticle.h"

int FlattenedVectorArray::rows() {
	return size() / cols();
}

int FlattenedVectorArray::cols() {
	return 3;
}

size_t FlattenedVectorArray::size() {
	return data.size();
}

void FlattenedVectorArray::set(int p_idx, LR_vector &v) {
	int base_idx = p_idx * cols();
	data[base_idx] = v.x;
	data[base_idx + 1] = v.y;
	data[base_idx + 2] = v.z;
}

void FlattenedConfigInfo::update(long long int step, const std::vector<BaseParticle *> &particles) {
	if(step == last_updated) {
		return;
	}

	last_updated = step;
	int N = particles.size();

	positions.data.resize(N * positions.cols());
	a1s.data.resize(N * a1s.cols());
	a3s.data.resize(N * a3s.cols());
	types.resize(N);

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];

		positions.set(i, p->pos);
		a1s.set(i, p->orientationT.v1);
		a3s.set(i, p->orientationT.v3);
		types[i] = p->type;
	}
}

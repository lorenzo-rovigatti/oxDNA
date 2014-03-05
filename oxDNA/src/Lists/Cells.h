/*
 * Cells.h
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#ifndef CELLS_H_
#define CELLS_H_

#include "BaseList.h"

/**
 * @brief Implementation of simple simulation cells.
 *
 * Internally, this class uses linked-list to keep track of particles in order to
 * have a computational complexity of O(N).
 */
template<typename number>
class Cells: public BaseList<number> {
protected:
	BaseParticle<number> **_heads;
	BaseParticle<number> **_next;
	int *_cells;
	int _N_cells;
	int _N_cells_side;
	number _sqr_rcut;

	inline int _get_cell_index(const LR_vector<number> &pos);
public:
	Cells(int &N, number &box);
	virtual ~Cells();

	virtual void init(BaseParticle<number> **particles, number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle<number> *p);
	virtual void global_update();
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p);
};

template<>
inline int Cells<float>::_get_cell_index(const LR_vector<float> &pos) {
	int res = (int) ((pos.x / this->_box - floorf(pos.x / this->_box)) * (1.f - FLT_EPSILON) * _N_cells_side);
	res += _N_cells_side * ((int) ((pos.y / this->_box - floorf(pos.y / this->_box)) * (1.f - FLT_EPSILON) * _N_cells_side));
	res += _N_cells_side * _N_cells_side * ((int) ((pos.z / this->_box - floorf(pos.z / this->_box)) * (1.f - FLT_EPSILON) * _N_cells_side));
	return res;
}

template<>
inline int Cells<double>::_get_cell_index(const LR_vector<double> &pos) {
	int res = (int) ((pos.x / this->_box - floor(pos.x / this->_box)) * (1. - DBL_EPSILON) * _N_cells_side);
	res += _N_cells_side * ((int) ((pos.y / this->_box - floor(pos.y / this->_box)) * (1. - DBL_EPSILON) * _N_cells_side));
	res += _N_cells_side * _N_cells_side * ((int) ((pos.z / this->_box - floor(pos.z / this->_box)) * (1. - DBL_EPSILON) * _N_cells_side));
	return res;
}

#endif /* CELLS_H_ */

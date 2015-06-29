/*
 * Cells.h
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#ifndef CELLS_H_
#define CELLS_H_

#include "BaseList.h"
#include <cfloat>

/**
 * @brief Implementation of simple simulation cells.
 *
 * Internally, this class uses linked-list to keep track of particles in order to
 * have a computational complexity of O(N).
 */
template<typename number>
class Cells: public BaseList<number> {
protected:
	int _allowed_type;
	bool _unlike_type_only;
	BaseParticle<number> **_heads;
	BaseParticle<number> **_next;
	int *_cells;
	int _N_cells;
	int _N_cells_side[3];
	number _sqr_rcut;
	LR_vector<number> _box_sides;

	inline int _get_cell_index(const LR_vector<number> &pos);
public:
	Cells(int &N, BaseBox<number> *box);
	virtual ~Cells();

	virtual void init(BaseParticle<number> **particles, number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle<number> *p);
	virtual void global_update(bool force_update=false);
	virtual std::vector<BaseParticle<number> *> get_neigh_list(BaseParticle<number> *p, bool all=false);
	virtual void change_box();

	virtual void set_allowed_type(int type) { _allowed_type = type; }
	virtual void set_unlike_type_only() { _unlike_type_only = true; }
};

template<>
inline int Cells<float>::_get_cell_index(const LR_vector<float> &pos) {
	int res = (int) ((pos.x / _box_sides.x - floorf(pos.x / _box_sides.x)) * (1.f - FLT_EPSILON) * _N_cells_side[0]);
	res += _N_cells_side[0] * ((int) ((pos.y / _box_sides.y - floorf(pos.y / _box_sides.y)) * (1.f - FLT_EPSILON) * _N_cells_side[1]));
	res += _N_cells_side[0] * _N_cells_side[1] * ((int) ((pos.z / _box_sides.z - floorf(pos.z / _box_sides.z)) * (1.f - FLT_EPSILON) * _N_cells_side[2]));
	return res;
}

template<>
inline int Cells<double>::_get_cell_index(const LR_vector<double> &pos) {
	int res = (int) ((pos.x / _box_sides[0] - floor(pos.x / _box_sides.x)) * (1. - DBL_EPSILON) * _N_cells_side[0]);
	res += _N_cells_side[0] * ((int) ((pos.y / _box_sides.y - floor(pos.y / _box_sides.y)) * (1. - DBL_EPSILON) * _N_cells_side[1]));
	res += _N_cells_side[0] * _N_cells_side[1] * ((int) ((pos.z / _box_sides.z - floor(pos.z / _box_sides.z)) * (1. - DBL_EPSILON) * _N_cells_side[2]));
	return res;
}

#endif /* CELLS_H_ */

/*
 * Cells.h
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#ifndef CELLS_H_
#define CELLS_H_

#include "BaseList.h"

#include <limits>

/**
 * @brief Implementation of simple simulation cells.
 *
 * Internally, this class uses linked-list to keep track of particles in order to
 * have a computational complexity of O(N).
 */

class Cells: public BaseList {
protected:
	int _allowed_type;
	bool _unlike_type_only;
	std::vector<BaseParticle *> _heads;
	std::vector<BaseParticle *> _next;
	int *_cells;
	int _N_cells;
	int _N_cells_side[3];
	number _sqr_rcut;
	bool _auto_optimisation;
	bool _lees_edwards;
	number _shear_rate;
	number _dt;

	void _set_N_cells_side_from_box(int N_cells_side[3], BaseBox *box);
	std::vector<BaseParticle *> _get_neigh_list(BaseParticle *p, bool all);
public:
	Cells(std::vector<BaseParticle *> &ps, BaseBox *box);
	Cells() = delete;
	virtual ~Cells();

	virtual void get_settings(input_file &inp);
	virtual void init(number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle *p);
	virtual void global_update(bool force_update=false);
	virtual std::vector<BaseParticle *> get_neigh_list(BaseParticle *p);
	virtual std::vector<BaseParticle *> get_complete_neigh_list(BaseParticle *p);

	virtual void set_allowed_type(int type) { _allowed_type = type; }
	virtual void set_unlike_type_only() { _unlike_type_only = true; }

	virtual int get_N_cells() { return _N_cells; }
	inline int get_cell_index(const LR_vector &pos);
};

inline int Cells::get_cell_index(const LR_vector &pos) {
	int res = (int) ((pos.x / this->_box_sides.x - floor(pos.x / this->_box_sides.x)) * (1. - std::numeric_limits<number>::epsilon())*_N_cells_side[0]);
	res += _N_cells_side[0]*((int) ((pos.y / this->_box_sides.y - floor(pos.y / this->_box_sides.y))*(1. - std::numeric_limits<number>::epsilon())*_N_cells_side[1]));
	res += _N_cells_side[0]*_N_cells_side[1]*((int) ((pos.z / this->_box_sides.z - floor(pos.z / this->_box_sides.z))*(1. - std::numeric_limits<number>::epsilon())*_N_cells_side[2]));
	return res;
}

#endif /* CELLS_H_ */

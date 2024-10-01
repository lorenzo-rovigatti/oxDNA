/*
 * RodCells.h
 *
 *  Created on: 7/oct/2018
 *      Author: flavio
 */

#ifndef RODCELLS_H_
#define RODCELLS_H_

#include "BaseList.h"

#include <limits>
#include <vector>
//#include <set>
//#include <unordered_set>
//#include <list>

/**
 * @brief Implementation of simple simulation cells.
 *
 * Internally, this class uses linked-list to keep track of particles in order to
 * have a computational complexity of O(N).
 */

class RodCells: public BaseList {
protected:
	int * _heads;
	int * _next;
	int *_cells;
	int _N_cells;
	int _N_cells_side[3];

	size_t _max_size;

	number _sqr_rcut;

	number _rod_cell_rcut;
	number _rod_length, _rod_length_2;
	int * _n_virtual_sites;
	int _n_part_types;
	int _n_virtual_sites_max;
	int _restrict_to_type;

	//std::set<int> _neighs;
	//std::unordered_set<int> _neighs; // unordered_set are measurably and consistently faster than the other options, but slower than sorting the resulting vector
	//std::list<int> _neighs;
	std::vector<bool> _added;

	void _set_N_cells_side_from_box(int N_cells_side[3], BaseBox *box);
	std::vector<BaseParticle *> _get_neigh_list(BaseParticle *p, bool all);
public:
	RodCells(std::vector<BaseParticle *> &ps, BaseBox *box);
	RodCells() = delete;
	virtual ~RodCells();

	virtual void get_settings(input_file &inp);
	virtual void init(number rcut);

	virtual bool is_updated();
	virtual void single_update(BaseParticle *p);
	virtual void global_update(bool force_update=false);
	virtual std::vector<BaseParticle *> get_neigh_list(BaseParticle *p);
	virtual std::vector<BaseParticle *> get_complete_neigh_list(BaseParticle *p);

	std::vector<BaseParticle * > whos_there(int idx);
	
	virtual int get_N_cells() { return _N_cells; }
	inline int get_cell_index(const LR_vector &pos);
};

inline int RodCells::get_cell_index(const LR_vector &pos) {
	int res = (int) ((pos.x / this->_box_sides.x - floor(pos.x / this->_box_sides.x)) * (1. - std::numeric_limits<number>::epsilon())*_N_cells_side[0]);
	res += _N_cells_side[0]*((int) ((pos.y / this->_box_sides.y - floor(pos.y / this->_box_sides.y))*(1. - std::numeric_limits<number>::epsilon())*_N_cells_side[1]));
	res += _N_cells_side[0]*_N_cells_side[1]*((int) ((pos.z / this->_box_sides.z - floor(pos.z / this->_box_sides.z))*(1. - std::numeric_limits<number>::epsilon())*_N_cells_side[2]));
	return res;
}

#endif /* RODCELLS_H_ */

/*
 * Cells.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "Cells.h"

template<typename number>
Cells<number>::Cells(int &N, number &box) : BaseList<number>(N, box) {
	_heads = NULL;
	_next = NULL;
	_cells = NULL;
	_N_cells = 0;
	_N_cells_side = 0;
	_sqr_rcut = 0;
}

template<typename number>
Cells<number>::~Cells() {
	if(_heads != NULL) delete[] _heads;
	if(_next != NULL) delete[] _next;
	if(_cells != NULL) delete[] _cells;
}

template<typename number>
void Cells<number>::init(BaseParticle<number> **particles, number rcut) {
	BaseList<number>::init(particles, rcut);

	_sqr_rcut = rcut*rcut;

	global_update();

	OX_LOG(Logger::LOG_INFO, "N_cells_side: %d", _N_cells_side);
}

template<typename number>
bool Cells<number>::is_updated() {
	int new_N_cells_side = (int) floor(this->_box / this->_rcut + 0.1);
	return (new_N_cells_side == _N_cells_side);
}

template<typename number>
void Cells<number>::single_update(BaseParticle<number> *p) {
	int old_cell = _cells[p->index];
	int new_cell = _get_cell_index(p->pos);

	if(old_cell != new_cell) {
		// remove p from its old cell
		BaseParticle<number> *previous = P_VIRTUAL;
		BaseParticle<number> *current = _heads[old_cell];
		while(current != p) {
			previous = current;
			current = _next[current->index];
		}
		if(previous == P_VIRTUAL) _heads[old_cell] = _next[p->index];
		else _next[previous->index] = _next[p->index];

		// add it to the new cell
		_next[p->index] = _heads[new_cell];
		_heads[new_cell] = p;

		_cells[p->index] = new_cell;
	}
}

template<typename number>
void Cells<number>::global_update() {
	_N_cells_side = (int) floor(this->_box / this->_rcut + 0.1);
	while(_N_cells_side > 750) _N_cells_side--;

	//if(_N_cells_side < 3) throw oxDNAException("%s, N_cells_side (%d) must be > 2", __FILE__, _N_cells_side);
	if (_N_cells_side < 3) _N_cells_side = 3;

	_N_cells = _N_cells_side * _N_cells_side * _N_cells_side;

	// this is needed only at the very beginning
	if(_heads != NULL) {
		delete[] _heads;
		delete[] _next;
		delete[] _cells;
	}
	_heads = new BaseParticle<number> *[_N_cells];
	_next = new BaseParticle<number> *[this->_N];
	_cells = new int[this->_N];

	for(int i = 0; i < _N_cells; i++) _heads[i] = P_VIRTUAL;
	for(int i = 0; i < this->_N; i++) _next[i] = P_VIRTUAL;

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		int cell_index = _get_cell_index(p->pos);
		BaseParticle<number> *old_head = _heads[cell_index];
		_heads[cell_index] = p;
		_cells[i] = cell_index;
		_next[i] = old_head;
	}
}

template<typename number>
std::vector<BaseParticle<number> *> Cells<number>::get_neigh_list(BaseParticle<number> *p) {
	std::vector<BaseParticle<number> *> res;

	int cind = _cells[p->index];
	int ind[3] = {
		cind % _N_cells_side,
		(cind / _N_cells_side) % _N_cells_side,
		cind / (_N_cells_side * _N_cells_side)
	};
	int loop_ind[3];
	for(int j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + _N_cells_side) % _N_cells_side;
		for(int k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + _N_cells_side) % _N_cells_side;
			for(int l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + _N_cells_side) % _N_cells_side;
				int loop_index = loop_ind[0] + _N_cells_side*(loop_ind[1] + loop_ind[2]*_N_cells_side);

				BaseParticle<number> *q = _heads[loop_index];
				while(q != P_VIRTUAL) {
					// if this is an MC simulation then we need full lists, otherwise the i-th particle will have neighbours with index > i
					if((p->index < q->index || (this->_is_MC && p != q)) && !p->is_bonded(q) && p->pos.sqr_min_image_distance(q->pos, this->_box) < _sqr_rcut) {
						res.push_back(q);
					}

					q = _next[q->index];
				}
			}
		}
	}

	return res;
}

template class Cells<float>;
template class Cells<double>;

/*
 * Cells.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "Cells.h"

#include "../Utilities/ConfigInfo.h"

Cells::Cells(std::vector<BaseParticle *> &ps, BaseBox *box) :
				BaseList(ps, box) {
	_cells = NULL;
	_N_cells = 0;
	_N_cells_side[0] = _N_cells_side[1] = _N_cells_side[2] = 0;
	_sqr_rcut = 0;
	_allowed_type = -1;
	_unlike_type_only = false;
	_auto_optimisation = true;
	_lees_edwards = false;
	_shear_rate = 0.;
	_dt = 0.;
}

Cells::~Cells() {
	if(_cells != NULL) {
		delete[] _cells;
	}
}

void Cells::get_settings(input_file &inp) {
	BaseList::get_settings(inp);
	getInputBool(&inp, "cells_auto_optimisation", &_auto_optimisation, 0);
	getInputBool(&inp, "lees_edwards", &_lees_edwards, 0);
	getInputNumber(&inp, "lees_edwards_shear_rate", &_shear_rate, 0);
	getInputNumber(&inp, "dt", &_dt, 0);
}

void Cells::init(number rcut) {
	BaseList::init(rcut);

	_sqr_rcut = rcut * rcut;

	global_update(true);

	OX_LOG(Logger::LOG_INFO, "(Cells.cpp) N_cells_side: %d, %d, %d; rcut=%g, IS_MC: %d", _N_cells_side[0], _N_cells_side[1], _N_cells_side[2], this->_rcut, this->_is_MC);
}

void Cells::_set_N_cells_side_from_box(int N_cells_side[3], BaseBox *box) {
	LR_vector box_sides = box->box_sides();
	number max_factor = pow(2. * _particles.size() / box->V(), 1. / 3.);
	for(int i = 0; i < 3; i++) {
		N_cells_side[i] = (int) (floor(box_sides[i] / this->_rcut) + 0.1);
		if(_auto_optimisation && N_cells_side[i] > ceil(max_factor * box_sides[i])) N_cells_side[i] = ceil(max_factor * box_sides[i]);

		if(N_cells_side[i] < 3) N_cells_side[i] = 3;
	}
}

bool Cells::is_updated() {
	int new_N_cells_side[3];
	_set_N_cells_side_from_box(new_N_cells_side, this->_box);
	return (new_N_cells_side[0] == _N_cells_side[0] && new_N_cells_side[1] == _N_cells_side[1] && new_N_cells_side[2] == _N_cells_side[2]);
}

void Cells::single_update(BaseParticle *p) {
	int old_cell = _cells[p->index];
	int new_cell = get_cell_index(p->pos);

	if(old_cell != new_cell) {
		// remove p from its old cell
		BaseParticle *previous = P_VIRTUAL;
		BaseParticle *current = _heads[old_cell];
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

void Cells::global_update(bool force_update) {
	this->_box_sides = this->_box->box_sides();
	_set_N_cells_side_from_box(_N_cells_side, this->_box);
	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	// we have to deallocate only if this is not the first time that this method gets called
	if(_heads.size() > 0) {
		_heads.clear();
		_next.clear();
		delete[] _cells;
	}
	_heads.resize(_N_cells, P_VIRTUAL);
	_next.resize(_particles.size(), P_VIRTUAL);
	_cells = new int[_particles.size()];

	for(uint i = 0; i < _particles.size(); i++) {
		BaseParticle *p = this->_particles[i];
		if(_allowed_type == -1 || p->type == _allowed_type) {
			int cell_index = get_cell_index(p->pos);
			if(cell_index < 0 || cell_index > _N_cells || std::isnan(cell_index) || std::isinf(cell_index)) {
				throw oxDNAException("Invalid cell %d for particle %d (pos: %lf %lf %lf)", cell_index, p->index, p->pos[0], p->pos[1], p->pos[2]);
			}
			BaseParticle *old_head = _heads[cell_index];
			_heads[cell_index] = p;
			_cells[i] = cell_index;
			_next[i] = old_head;
		}
	}
}

std::vector<BaseParticle *> Cells::_get_neigh_list(BaseParticle *p, bool all) {
	// the static here makes sure that the memory allocated every time the function is called is not free'd
	// so that only push_back's that make the vector go beyond their last size trigger re-allocation
	static std::vector<BaseParticle *> res;
	res.clear();

	int cind = _cells[p->index];
	int ind[3] = { cind % _N_cells_side[0], (cind / _N_cells_side[0]) % _N_cells_side[1], cind / (_N_cells_side[0] * _N_cells_side[1]) };
	int loop_ind[3];

	// y direction
	for(int k = -1; k < 2; k++) {
		loop_ind[1] = (ind[1] + k + _N_cells_side[1]) % _N_cells_side[1];

		int cx_lower = -1;
		int cx_upper = 2;
		if(_lees_edwards && _N_cells_side[0] > 3) {
			bool lower_edge = (ind[1] == 0) && k == -1;
			bool upper_edge = (ind[1] == (_N_cells_side[1] - 1)) && k == 1;
			if(lower_edge || upper_edge) {
				const LR_vector &L = this->_box->box_sides();
				number delta_x = _shear_rate * L.y * _dt * CONFIG_INFO->curr_step;
				delta_x -= floor(delta_x / L.x) * L.x;
				int c_delta_x = (delta_x / L.x) * _N_cells_side[0];
				if(lower_edge) {
					cx_lower += c_delta_x;
					cx_upper = cx_lower + 4;
				}
				if(upper_edge) {
					cx_lower -= c_delta_x + 1;
					cx_upper = cx_lower + 4;
				}
			}
		}

		// x direction
		for(int j = cx_lower; j < cx_upper; j++) {
			// the factor 2 in the bracketed expression is required when Lees-Edwards conditions are enabled
			loop_ind[0] = (ind[0] + j + 2 * _N_cells_side[0]) % _N_cells_side[0];

			// z direction
			for(int l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + _N_cells_side[2]) % _N_cells_side[2];
				int loop_index = loop_ind[0] + _N_cells_side[0] * (loop_ind[1] + loop_ind[2] * _N_cells_side[1]);

				BaseParticle *q = _heads[loop_index];
				while(q != P_VIRTUAL) {
					// if this is an MC simulation or all == true we need full lists, otherwise the i-th particle will have neighbours with index > i
					bool include_q = (p != q) && (all || ((p->index > q->index || this->_is_MC)));
					include_q = include_q && (!_unlike_type_only || p->type != q->type);
					if(include_q && !p->is_bonded(q) && this->_box->sqr_min_image_distance(p->pos, q->pos) < _sqr_rcut) {
						res.push_back(q);
					}

					q = _next[q->index];
				}
			}
		}
	}

	return res;
}

std::vector<BaseParticle *> Cells::get_neigh_list(BaseParticle *p) {
	return _get_neigh_list(p, false);
}

std::vector<BaseParticle *> Cells::get_complete_neigh_list(BaseParticle *p) {
	return _get_neigh_list(p, true);
}

/*
 * Cells.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "Cells.h"

#include "../Utilities/ConfigInfo.h"

template<typename number>
Cells<number>::Cells(int &N, BaseBox<number> *box) : BaseList<number>(N, box) {
	_heads = NULL;
	_next = NULL;
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

template<typename number>
Cells<number>::~Cells() {
	if(_heads != NULL) delete[] _heads;
	if(_next != NULL) delete[] _next;
	if(_cells != NULL) delete[] _cells;
}

template<typename number>
void Cells<number>::get_settings(input_file &inp) {
	BaseList<number>::get_settings(inp);
	getInputBool(&inp, "cells_auto_optimisation", &_auto_optimisation, 0);
	getInputBool(&inp, "lees_edwards", &_lees_edwards, 0);
	getInputNumber(&inp, "lees_edwards_shear_rate", &_shear_rate, 0);
	getInputNumber(&inp, "dt", &_dt, 0);
}

template<typename number>
void Cells<number>::init(BaseParticle<number> **particles, number rcut) {
	BaseList<number>::init(particles, rcut);

	_sqr_rcut = rcut*rcut;

	global_update(true);

	OX_LOG(Logger::LOG_INFO, "(Cells.cpp) N_cells_side: %d, %d, %d; rcut=%g, IS_MC: %d", _N_cells_side[0], _N_cells_side[1], _N_cells_side[2], this->_rcut, this->_is_MC);
}

/*
template<typename number>
void Cells<number>::dump(std::string filename) {
	std::ofstream output;
	output.open(filename.c_str(), std::ifstream::out | std::ifstream::binary);

	if (!output.good()) throw oxDNAException("Could not open file %s to which to write lists.", filename.c_str());

	int minus1 = -1;

	output.write ((char *) &this->_N, sizeof(int));
	output.write ((char *) _N_cells_side, 3 * sizeof(int));
	
	for (int icell = 0; icell < _N_cells; icell ++) {
		BaseParticle<number> * p = _heads[icell];
		if (p != P_VIRTUAL) output.write ((char *) &(p->index), sizeof(int));
		else output.write ((char *) &minus1, sizeof(int));
	}

	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> * p = _next[i];
		if (p != P_VIRTUAL) output.write ((char *) &(p->index), sizeof(int));
		else output.write ((char *) &minus1, sizeof(int));
	}
	
	output.close();
}

template<typename number>
void Cells<number>::load(std::string filename) {
	std::ifstream input;
	input.open(filename.c_str(), std::ifstream::in | std::ifstream::binary);

	if (!input.good()) throw oxDNAException("Could not open file %s from which to load lists.", filename.c_str());

	_box_sides = this->_box->box_sides();

	int tmpi;

	// we have to get the number of particles
	input.read((char *)(&(this->_N)), sizeof(int));
	
	// we have to get the number of cells
	input.read((char *) _N_cells_side, 3 * sizeof(int));
	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	_heads = new BaseParticle<number> *[_N_cells];
	_next = new BaseParticle<number> *[this->_N];
	_cells = new int[this->_N];

	// read the cells' heads
	for(int i = 0; i < _N_cells; i ++) {
		input.read((char *) &tmpi, sizeof(int));
		if (tmpi > this->_N) throw oxDNAException("Error reading lists: found index %d, which is larger than N=%d", tmpi, this->_N);
		if (tmpi >= 0) _heads[i] = this->_particles[tmpi];
		else _heads[i] = P_VIRTUAL;
	}

	// read next for each particle
	for(int i = 0; i < this->_N; i ++) {
		input.read((char *) &tmpi, sizeof(int));
		if (tmpi > this->_N) throw oxDNAException("Error reading lists: found index %d, which is larger than N=%d", tmpi, this->_N);
		if (tmpi >= 0) _next[i] = this->_particles[tmpi];
		else _next[i] = P_VIRTUAL;
	}
	
	// now we just need to set the cell index
	for(int i = 0; i < _N_cells; i ++) {
		BaseParticle<number> * p = _heads[i];
		while (p != P_VIRTUAL) {
			_cells[p->index] = i;
			p = _next[p->index];
		}
	}

	input.close();
}
*/

template<typename number>
void Cells<number>::_set_N_cells_side_from_box(int N_cells_side[3], BaseBox<number> *box) {
	LR_vector<number> box_sides = box->box_sides();
	number max_factor = pow(2.*this->_N/box->V(), 1./3.);
	for(int i = 0; i < 3; i++) {
		N_cells_side[i] = (int) (floor(box_sides[i] / this->_rcut) + 0.1);
		if(_auto_optimisation && N_cells_side[i] > ceil(max_factor*box_sides[i])) N_cells_side[i] = ceil(max_factor*box_sides[i]);

		if(N_cells_side[i] < 3) N_cells_side[i] = 3;
	}
}

template<typename number>
bool Cells<number>::is_updated() {
	int new_N_cells_side[3];
	_set_N_cells_side_from_box(new_N_cells_side, this->_box);
	return (new_N_cells_side[0] == _N_cells_side[0] && new_N_cells_side[1] == _N_cells_side[1] && new_N_cells_side[2] == _N_cells_side[2]);
}

template<typename number>
void Cells<number>::single_update(BaseParticle<number> *p) {
	int old_cell = _cells[p->index];
	int new_cell = get_cell_index(p->pos);

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
void Cells<number>::global_update(bool force_update) {
	this->_box_sides = this->_box->box_sides();
	_set_N_cells_side_from_box(_N_cells_side, this->_box);
	_N_cells = _N_cells_side[0]*_N_cells_side[1]*_N_cells_side[2];

	// we have to deallocate only if this is not the first time that this method gets called
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
		if(_allowed_type == -1 || p->type == _allowed_type) {
			int cell_index = get_cell_index(p->pos);
			BaseParticle<number> *old_head = _heads[cell_index];
			_heads[cell_index] = p;
			_cells[i] = cell_index;
			_next[i] = old_head;
		}
	}
}

template<typename number>
std::vector<BaseParticle<number> *> Cells<number>::_get_neigh_list(BaseParticle<number> *p, bool all) {
	std::vector<BaseParticle<number> *> res;

	int cind = _cells[p->index];
	int ind[3] = {
		cind % _N_cells_side[0],
		(cind / _N_cells_side[0]) % _N_cells_side[1],
		cind / (_N_cells_side[0]*_N_cells_side[1])
	};
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
				const LR_vector<number> &L = this->_box->box_sides();
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
			loop_ind[0] = (ind[0] + j + 2*_N_cells_side[0]) % _N_cells_side[0];

			// z direction
			for(int l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + _N_cells_side[2]) % _N_cells_side[2];
				int loop_index = loop_ind[0] + _N_cells_side[0]*(loop_ind[1] + loop_ind[2]*_N_cells_side[1]);

				BaseParticle<number> *q = _heads[loop_index];
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

template<typename number>
std::vector<BaseParticle<number> *> Cells<number>::get_neigh_list(BaseParticle<number> *p) {
	return _get_neigh_list(p, false);
}

template<typename number>
std::vector<BaseParticle<number> *> Cells<number>::get_complete_neigh_list(BaseParticle<number> *p) {
	return _get_neigh_list(p, true);
}

template class Cells<float>;
template class Cells<double>;

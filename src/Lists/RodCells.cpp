/*
 * RodCells.cpp
 *
 *  Created on: 05/oct/2018
 *      Author: flavio
 */

#include "RodCells.h"

#include "../Utilities/ConfigInfo.h"
#include <algorithm>

template<typename number>
RodCells<number>::RodCells(int &N, BaseBox<number> *box) : BaseList<number>(N, box) {
	_heads = NULL;
	_next = NULL;
	_cells = NULL;
	_N_cells = 0;
	_N_cells_side[0] = _N_cells_side[1] = _N_cells_side[2] = 0;
	_sqr_rcut = 0;
	_rod_cell_rcut = (number) 0.f;
	_rod_length = (number) -1.f;
	_n_virtual_sites = -1;
	_max_size = 20;
	_added = std::vector<bool> ();
	//_neighs = std::unordered_set<int> ();
	//_neighs = std::set<int> ();
	//_neighs = std::list<int> ();
}

template<typename number>
RodCells<number>::~RodCells() {
	if(_heads != NULL) delete[] _heads;
	if(_next != NULL) delete[] _next;
	if(_cells != NULL) delete[] _cells;
}

template<typename number>
void RodCells<number>::get_settings(input_file &inp) {
	BaseList<number>::get_settings(inp);
	getInputNumber(&inp, "rod_cell_rcut", &_rod_cell_rcut, 1);
	getInputNumber(&inp, "rod_length", &_rod_length, 1);
}

template<typename number>
void RodCells<number>::init(BaseParticle<number> **particles, number rcut) {
	BaseList<number>::init(particles, rcut);

	_sqr_rcut = rcut * rcut;

	number du = (2./sqrt(3.)) * _rod_cell_rcut;

	_n_virtual_sites = (int)(floor(_rod_length / du)) + 1;

	//printf ("@@ du, rcut: %g %g\n", du, _rod_cell_rcut);

	for (int i = 0; i < this->_N; i ++) _added.push_back(false);

	global_update(true);

	OX_LOG(Logger::LOG_INFO, "(RodCells.cpp) N_cells_side: %d, %d, %d; rcut=%g, rod_cell_size=%g, _n_virtual_sites=%d, IS_MC: %d",_N_cells_side[0], _N_cells_side[1], _N_cells_side[2], this->_rcut, _rod_cell_rcut, _n_virtual_sites, this->_is_MC);
}

template<typename number>
void RodCells<number>::_set_N_cells_side_from_box(int N_cells_side[3], BaseBox<number> *box) {
	LR_vector<number> box_sides = box->box_sides();  // TODO: perhaps use pointer instead of copying?
	for(int i = 0; i < 3; i++) {
		N_cells_side[i] = (int) (floor(box_sides[i] / _rod_cell_rcut) + 0.1);
		if(N_cells_side[i] < 3)
			N_cells_side[i] = 3;
	}
}

template<typename number>
bool RodCells<number>::is_updated() {
	int new_N_cells_side[3];
	_set_N_cells_side_from_box(new_N_cells_side, this->_box);
	return (new_N_cells_side[0] == _N_cells_side[0] && new_N_cells_side[1] == _N_cells_side[1] && new_N_cells_side[2] == _N_cells_side[2]);
}

template<typename number>
void RodCells<number>::single_update(BaseParticle<number> *p) {

	// remove from old cells
	for (int k = 0; k < _n_virtual_sites; k++) {
		int site_idx = p->index * _n_virtual_sites + k;
		int old_cell = _cells[site_idx];

		int previous = -1;
		int current = _heads[old_cell];
		while (current != site_idx) {
			previous = current;
			current = _next[current];
		}
		if (previous == -1) _heads[old_cell] = _next[site_idx];
		else _next[previous] = _next[site_idx];
	}

	// compute new cell indexes and add to new cells
	LR_vector<number> * r = &p->pos;
	LR_vector<number> * u = &p->orientation.v3;

	LR_vector<number> stride = (_rod_length / (_n_virtual_sites - 1)) * (*u);
	LR_vector<number> site_pos = (*r) - (_rod_length / (number) 2.f) * (*u);
	for (int k = 0; k < _n_virtual_sites; k ++) {
		int site_idx = p->index * _n_virtual_sites + k;
		int new_cell = get_cell_index(site_pos);
		_cells[site_idx] = new_cell;
		_next[site_idx] = _heads[new_cell];
		_heads[new_cell] = site_idx;

		site_pos += stride;
	}

	return;
}

template<typename number>
void RodCells<number>::global_update(bool force_update) {
	this->_box_sides = this->_box->box_sides();
	_set_N_cells_side_from_box(_N_cells_side, this->_box);
	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	// we have to deallocate only if this is not the first time that this method gets called
	if(_heads != NULL) {
		delete[] _heads;
		delete[] _next;
		delete[] _cells;
	}
	_heads = new int [_N_cells];
	_cells = new int [this->_N * _n_virtual_sites];
	_next = new int [this->_N * _n_virtual_sites];

	for(int i = 0; i < _N_cells; i++) _heads[i] = -1;
	for(int i = 0; i < this->_N * _n_virtual_sites; i++) _next[i] = -1;

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		LR_vector<number> stride = (_rod_length / (_n_virtual_sites - 1)) * p->orientation.v3;
		LR_vector<number> site_pos = p->pos - (_rod_length / (number) 2.f) * p->orientation.v3;
		for (int k = 0; k < _n_virtual_sites; k ++) {
			int site_idx = i * _n_virtual_sites + k;
			int cell_index = get_cell_index(site_pos);
			_cells[site_idx] = cell_index;
			_next[site_idx] = _heads[cell_index];
			_heads[cell_index] = site_idx;

			site_pos += stride;
		}
	}
}

template<typename number>
std::vector<BaseParticle<number> *> RodCells<number>::_get_neigh_list(BaseParticle<number> *p, bool all) {
	std::vector<BaseParticle<number> *> res;

	res.reserve(_max_size);

	//_neighs.clear();
	//typename std::unordered_set<int>::iterator it = _neighs.begin();
	//typename std::set<int>::iterator it = _neighs.begin();
	//typename std::vector<BaseParticle<number> *>::iterator it = res.begin();
	//typename std::list<int>::iterator it;

	int last_inserted = -4;

	for (int s = 0; s < _n_virtual_sites; s ++) {
		int site_idx = p->index * _n_virtual_sites + s;
		int cind = _cells[site_idx];
		int ind[3] = {
			cind % _N_cells_side[0],
			(cind / _N_cells_side[0]) % _N_cells_side[1],
			cind / (_N_cells_side[0]*_N_cells_side[1])
		};
		int loop_ind[3];

		for (int i = -1; i < 2; i ++) {
			loop_ind[0] = (ind[0] + i) % _N_cells_side[0];
			if (loop_ind[0] < 0) loop_ind[0] += _N_cells_side[0];
			for (int j = -1; j < 2; j ++) {
				loop_ind[1] = (ind[1] + j) % _N_cells_side[1];
				if (loop_ind[1] < 0) loop_ind[1] += _N_cells_side[1];
				for (int k = -1; k < 2; k ++) {
					loop_ind[2] = (ind[2] + k) % _N_cells_side[2];
					if (loop_ind[2] < 0) loop_ind[2] += _N_cells_side[2];
					int other_cell_index = loop_ind[0] + _N_cells_side[0]*(loop_ind[1] + _N_cells_side[1]*loop_ind[2]);

					int n = _heads[other_cell_index];

					/*
					// this cycle works well with sets
					while (n != -1) {
						//_neighs.insert(n / _n_virtual_sites); // simple version
						it = _neighs.insert(it, n / _n_virtual_sites);
						n = _next[n];
					}

					// cycle to work with lists
					while (n != -1) {
						_neighs.push_back(n / _n_virtual_sites);
						n = _next[n];
					}*/

					while (n != -1) {
						int m = n / _n_virtual_sites;
						if (m != p->index && m != last_inserted && _added[m] == false) {
							if (all || this->_is_MC || p->index > m) {
								res.push_back(this->_particles[m]);
								last_inserted = m;
								_added[m] = true;
							}
						}
						n = _next[n];
					}
				}
			}
		}
	}

	if (res.size() > _max_size) {
		// this will adjust the size of the reserved vector
		_max_size = res.size();
		//printf ("## new size: %d\n", (int) _max_size );
	}

	// sort the vector and then let unique do the magic
	//std::sort(res.begin(), res.end());
	//res.erase(std::unique(res.begin(), res.end()), res.end());

	for (auto it = res.begin(); it != res.end(); ++it) _added[(*it)->index] = false;

	/*
	 // this was needed for ordered sets
	for (it = _neighs.begin(); it != _neighs.end(); it ++) {
		if (*it == p->index)
			continue;
		if (all || this->_is_MC || p->index > *it) {
			//if (this->_box->sqr_min_image_distance(p->pos, this->_particles[*it]->pos) < _sqr_rcut)
			res.push_back(this->_particles[*it]);
		}
	}

	//lists
	_neighs.sort();
	last_inserted = -4;
	for (it = _neighs.begin(); it != _neighs.end(); it ++) {
		if (*it == p->index || *it == last_inserted)
			continue;

		if (all || this->_is_MC || p->index > *it) {
			//if (this->_box->sqr_min_image_distance(p->pos, this->_particles[*it]->pos) < _sqr_rcut)
			res.push_back(this->_particles[*it]);
			last_inserted = *it;
		}
	}
	*/

	return res;
}

template<typename number>
std::vector<BaseParticle<number> *> RodCells<number>::get_neigh_list(BaseParticle<number> *p) {
	return _get_neigh_list(p, false);
}

template<typename number>
std::vector<BaseParticle<number> *> RodCells<number>::get_complete_neigh_list(BaseParticle<number> *p) {
	return _get_neigh_list(p, true);
}

template class RodCells<float>;
template class RodCells<double>;

/*
 * RodCells.cpp
 *
 *  Created on: 05/oct/2018
 *      Author: flavio
 */

#include "RodCells.h"

#include "../Utilities/ConfigInfo.h"
#include <algorithm>

RodCells::RodCells(std::vector<BaseParticle *> &ps, BaseBox *box) :
				BaseList(ps, box) {
	_heads = NULL;
	_next = NULL;
	_cells = NULL;
	_N_cells = 0;
	_N_cells_side[0] = _N_cells_side[1] = _N_cells_side[2] = 0;
	_sqr_rcut = 0;
	_rod_cell_rcut = (number) 0.f;
	_rod_length = (number) -1.f;
	_rod_length_2 = (number) -1.f;
	_n_virtual_sites = NULL;
	_max_size = 20;
	_added = std::vector<bool>();
	_n_part_types = 1;
	_n_virtual_sites_max = -1;
	_restrict_to_type = -1;
	//_neighs = std::unordered_set<int> ();
	//_neighs = std::set<int> ();
	//_neighs = std::list<int> ();
}

RodCells::~RodCells() {
	if(_heads != NULL) delete[] _heads;
	if(_next != NULL) delete[] _next;
	if(_cells != NULL) delete[] _cells;
	if(_n_virtual_sites != NULL) delete[] _n_virtual_sites;
}

void RodCells::get_settings(input_file &inp){
	BaseList::get_settings(inp);
	getInputNumber(&inp, "rod_cell_rcut", &_rod_cell_rcut, 1);
	getInputNumber(&inp, "rod_length", &_rod_length, 1);
	getInputNumber(&inp, "rod_length_2", &_rod_length_2, 1);
	getInputInt(&inp, "rod_cell_n_part_types", &_n_part_types, 0);
	getInputInt(&inp, "rod_cell_restrict_to_type", &_restrict_to_type, 0);
}

void RodCells::init(number rcut) {
	BaseList::init(rcut);

	_sqr_rcut = rcut * rcut;

	_n_virtual_sites = new int[_n_part_types];
	_n_virtual_sites[0] = 1;
	while((_rod_length / _n_virtual_sites[0]) + (number) 0.01 >= _rod_cell_rcut)
		_n_virtual_sites[0] += 1;

	if(_n_part_types == 2) {
		if (_rod_length_2 < 0.) throw oxDNAException("(RodCells.cpp) Can't run with 2 particle types and no _rod_length_2");
		_n_virtual_sites[1] = 1;
		while((_rod_length_2 / _n_virtual_sites[1]) + (number) 0.01 >= _rod_cell_rcut)
			_n_virtual_sites[1] += 1;
	}

	if(_n_part_types > 2) throw oxDNAException("RodCells.ccp can handle at most 2 particle types for now");

	for(uint i = 0; i < _particles.size(); i++) {
		if(_particles[i]->type >= _n_part_types) throw oxDNAException("Found particle with index %d and type %d, but RodCell is set up with only %d particle types", _particles[i]->index, _particles[i]->type, _n_part_types);
	}

	_n_virtual_sites_max = _n_virtual_sites[0] > _n_virtual_sites[1] ? _n_virtual_sites[0] : _n_virtual_sites[1];

	for(uint i = 0; i < _particles.size(); i++)
		_added.push_back(false);

	global_update(true);

	OX_LOG(Logger::LOG_INFO, "(RodCells.cpp) N_cells_side: %d, %d, %d; rcut=%g, rod_cell_size=%g, _n_virtual_sites=%d, restrict_to_type=%d, IS_MC: %d, stride: %g",_N_cells_side[0], _N_cells_side[1], _N_cells_side[2], this->_rcut, _rod_cell_rcut, _n_virtual_sites[0], _restrict_to_type, this->_is_MC, _rod_length/_n_virtual_sites[0]);
}

void RodCells::_set_N_cells_side_from_box(int N_cells_side[3], BaseBox *box) {
	LR_vector box_sides = box->box_sides();  // TODO: perhaps use pointer instead of copying?
	number my_rcut = _rod_cell_rcut + _rod_length / 2.f / _n_virtual_sites[0] + _rod_length_2 / 2.f / _n_virtual_sites[1];
	//OX_LOG(Logger::LOG_INFO, "(RodCells.cpp) myrcut: %g (%g %g)", my_rcut, + _rod_length / 2.f / _n_virtual_sites[0], + _rod_length_2 / 2.f / _n_virtual_sites[1]);
	for(int i = 0; i < 3; i++) {
		N_cells_side[i] = (int) (floor(box_sides[i] / my_rcut) + 0.1);
		if(N_cells_side[i] < 3) N_cells_side[i] = 3;
	}
	while((N_cells_side[0] * N_cells_side[1] * N_cells_side[2]) > (int) (2 * _particles.size() * _n_virtual_sites_max)) {
		N_cells_side[0]--;
		N_cells_side[1]--;
		N_cells_side[2]--;
	}
}

bool RodCells::is_updated() {
	int new_N_cells_side[3];
	_set_N_cells_side_from_box(new_N_cells_side, this->_box);
	return (new_N_cells_side[0] == _N_cells_side[0] && new_N_cells_side[1] == _N_cells_side[1] && new_N_cells_side[2] == _N_cells_side[2]);
}

void RodCells::single_update(BaseParticle *p) {

	//if (_restrict_to_type >= 0 && p->type != _restrict_to_type) return;

	// remove from old cells
	for(int k = 0; k < _n_virtual_sites[p->type]; k++) {
		int site_idx = p->index * _n_virtual_sites_max + k;
		int old_cell = _cells[site_idx];

		int previous = -1;
		int current = _heads[old_cell];
		while(current != site_idx) {
			previous = current;
			current = _next[current];
		}
		if(previous == -1) _heads[old_cell] = _next[site_idx];
		else _next[previous] = _next[site_idx];
	}

	// compute new cell indexes and add to new cells
	LR_vector * r = &p->pos;
	LR_vector * u = &p->orientation.v3;
	number s = (_rod_length / (_n_virtual_sites[p->type]));
	LR_vector stride = s * (*u);
	LR_vector site_pos = (*r) - (0.5f * _rod_length + 0.5f * s) * (*u);
	if(_n_virtual_sites[p->type] < 2) site_pos = *r;
	for(int k = 0; k < _n_virtual_sites[p->type]; k++) {
		int site_idx = p->index * _n_virtual_sites_max + k;
		int new_cell = get_cell_index(site_pos);
		_cells[site_idx] = new_cell;
		_next[site_idx] = _heads[new_cell];
		_heads[new_cell] = site_idx;

		site_pos += stride;
	}

	return;
}

void RodCells::global_update(bool force_update) {
	this->_box_sides = this->_box->box_sides();
	_set_N_cells_side_from_box(_N_cells_side, this->_box);
	_N_cells = _N_cells_side[0] * _N_cells_side[1] * _N_cells_side[2];

	// we have to deallocate only if this is not the first time that this method gets called
	if(_heads != NULL) {
		delete[] _heads;
		delete[] _next;
		delete[] _cells;
	}
	_heads = new int[_N_cells];
	_cells = new int[_particles.size() * _n_virtual_sites_max];
	_next = new int[_particles.size() * _n_virtual_sites_max];

	for(int i = 0; i < _N_cells; i++) {
		_heads[i] = -1;
	}

	for(uint i = 0; i < _particles.size() * _n_virtual_sites_max; i++) {
		_next[i] = -1;
	}

	for(uint i = 0; i < _particles.size(); i++) {
		BaseParticle *p = this->_particles[i];
		//if (_restrict_to_type >= 0 && p->type != _restrict_to_type) continue;
		// FIX HERE with new stride
		LR_vector * r = &p->pos;
		LR_vector * u = &p->orientation.v3;
		number s = (_rod_length / (_n_virtual_sites[p->type]));
		LR_vector stride = s * (*u);
		LR_vector site_pos = (*r) - (0.5f * _rod_length + 0.5f * s) * (*u);
		if(_n_virtual_sites[p->type] < 2) site_pos = p->pos;
		for(int k = 0; k < _n_virtual_sites[p->type]; k++) {
			int site_idx = i * _n_virtual_sites_max + k;
			int cell_index = get_cell_index(site_pos);
			_cells[site_idx] = cell_index;
			_next[site_idx] = _heads[cell_index];
			_heads[cell_index] = site_idx;
			site_pos += stride;
		}
	}
}

std::vector<BaseParticle *> RodCells::_get_neigh_list(BaseParticle *p, bool all) {
	std::vector<BaseParticle *> res;

	int want_type = -1;
	if(_restrict_to_type >= 0) {
		if(p->type == _restrict_to_type) want_type = -1; // we want 'em all
		else want_type = _restrict_to_type;               // we just want the other type
	}

	res.reserve(_max_size);

	//_neighs.clear();
	//typename std::unordered_set<int>::iterator it = _neighs.begin();
	//typename std::set<int>::iterator it = _neighs.begin();
	//typename std::vector<BaseParticle *>::iterator it = res.begin();
	//typename std::list<int>::iterator it;

	int last_inserted = -4;

	for(int s = 0; s < _n_virtual_sites[p->type]; s++) {
		int site_idx = p->index * _n_virtual_sites_max + s;
		int cind = _cells[site_idx];
		int ind[3] = { cind % _N_cells_side[0], (cind / _N_cells_side[0]) % _N_cells_side[1], cind / (_N_cells_side[0] * _N_cells_side[1]) };
		int loop_ind[3];

		for(int i = -1; i < 2; i++) {
			loop_ind[0] = (ind[0] + i) % _N_cells_side[0];
			if(loop_ind[0] < 0) loop_ind[0] += _N_cells_side[0];
			for(int j = -1; j < 2; j++) {
				loop_ind[1] = (ind[1] + j) % _N_cells_side[1];
				if(loop_ind[1] < 0) loop_ind[1] += _N_cells_side[1];
				for(int k = -1; k < 2; k++) {
					loop_ind[2] = (ind[2] + k) % _N_cells_side[2];
					if(loop_ind[2] < 0) loop_ind[2] += _N_cells_side[2];
					int other_cell_index = loop_ind[0] + _N_cells_side[0] * (loop_ind[1] + _N_cells_side[1] * loop_ind[2]);

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

					while(n != -1) {
						int m = n / _n_virtual_sites_max;
						if(m != p->index && m != last_inserted) {
							if(all || this->_is_MC || p->index > m) {
								if(want_type < 0 || this->_particles[m]->type == want_type) {
									if(_added[m] == false) {
										res.push_back(this->_particles[m]);
										last_inserted = m;
										_added[m] = true;
									}
								}
							}
						}
						n = _next[n];
					}
				}
			}
		}
	}

	if(res.size() > _max_size) {
		// this will adjust the size of the reserved vector
		_max_size = res.size();
		//printf ("## new size: %d\n", (int) _max_size );
	}

	// sort the vector and then let unique do the magic
	//std::sort(res.begin(), res.end());
	//res.erase(std::unique(res.begin(), res.end()), res.end());

	typename std::vector<BaseParticle *>::iterator it;
	for(it = res.begin(); it != res.end(); ++it)
		_added[(*it)->index] = false;

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

std::vector<BaseParticle *> RodCells::whos_there(int idx) {
	if(idx >= _N_cells) throw oxDNAException("wrong cell idx");

	std::vector<BaseParticle *> res;
	res.reserve(_max_size);

	int ind[3] = { idx % _N_cells_side[0], (idx / _N_cells_side[0]) % _N_cells_side[1], idx / (_N_cells_side[0] * _N_cells_side[1]) };

	int loop_ind[3];
	for(int i = -1; i < 2; i++) {
		loop_ind[0] = (ind[0] + i) % _N_cells_side[0];
		if(loop_ind[0] < 0) loop_ind[0] += _N_cells_side[0];
		for(int j = -1; j < 2; j++) {
			loop_ind[1] = (ind[1] + j) % _N_cells_side[1];
			if(loop_ind[1] < 0) loop_ind[1] += _N_cells_side[1];
			for(int k = -1; k < 2; k++) {
				loop_ind[2] = (ind[2] + k) % _N_cells_side[2];
				if(loop_ind[2] < 0) loop_ind[2] += _N_cells_side[2];
				int other_cell_index = loop_ind[0] + _N_cells_side[0] * (loop_ind[1] + _N_cells_side[1] * loop_ind[2]);

				int n = _heads[other_cell_index];

				while(n != -1) {
					int m = n / _n_virtual_sites_max;
					if(_added[m] == false) {
						res.push_back(this->_particles[m]);
						_added[m] = true;
					}
					n = _next[n];
				}
			}
		}
	}

	typename std::vector<BaseParticle *>::iterator it;
	for(it = res.begin(); it != res.end(); ++it)
		_added[(*it)->index] = false;

	return res;
}

std::vector<BaseParticle *> RodCells::get_neigh_list(BaseParticle *p) {
	return _get_neigh_list(p, false);
}

std::vector<BaseParticle *> RodCells::get_complete_neigh_list(BaseParticle *p) {
	return _get_neigh_list(p, true);
}

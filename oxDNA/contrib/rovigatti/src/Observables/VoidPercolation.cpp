/*
 * VoidPercolation.cpp
 *
 *  Created on: Sep 27, 2016
 *      Author: Lorenzo Rovigatti
 */

#include "VoidPercolation.h"

#include "Boxes/BoxFactory.h"
#include "Utilities/Utils.h"

#include <sstream>

using namespace std;

template<typename number>
VPCells<number>::VPCells(int &N, BaseBox<number> *box) : Cells<number>(N, box) {
	_clusters.resize(this->_N, 0);
	_sizes.resize(this->_N + 1, 0);
	csd.resize(this->_N + 1, 0);
}

template<typename number>
VPCells<number>::~VPCells() {

}

template<typename number>
void VPCells<number>::_flip_neighs(BaseParticle<number> *p) {
	int p_c = _clusters[p->index];

	int cind = this->_cells[p->index];
	int ind[3] = {
		cind % this->_N_cells_side[0],
		(cind / this->_N_cells_side[0]) % this->_N_cells_side[1],
		cind / (this->_N_cells_side[0]*this->_N_cells_side[1])
	};
	int loop_ind[3];
	for(int j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + this->_N_cells_side[0]) % this->_N_cells_side[0];
		for(int k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + this->_N_cells_side[1]) % this->_N_cells_side[1];
			for(int l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + this->_N_cells_side[2]) % this->_N_cells_side[2];
				int loop_index = loop_ind[0] + this->_N_cells_side[0]*(loop_ind[1] + loop_ind[2]*this->_N_cells_side[1]);

				BaseParticle<number> *q = this->_heads[loop_index];
				while(q != P_VIRTUAL) {
					if(this->_box->sqr_min_image_distance(p->pos, q->pos) < this->_sqr_rcut) {
						int q_c = _clusters[q->index];
						if(q_c > p_c) {
							_clusters[q->index] = p_c;
							_flip_neighs(q);
						}
					}

					q = this->_next[q->index];
				}
			}
		}
	}
}

template<typename number>
void VPCells<number>::compute_csd() {
	fill(_clusters.begin(), _clusters.end(), 0);
	for(int i = 0; i < this->_N; i++) _clusters[i] = i;

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		_flip_neighs(p);
	}

	fill(_sizes.begin(), _sizes.end(), 0);
	for(vector<int>::iterator it = _clusters.begin(); it != _clusters.end(); it++) {
		int id = *it;
		_sizes[id]++;
	}

	fill(csd.begin(), csd.end(), 0);
	for(vector<int>::iterator it = _sizes.begin(); it != _sizes.end(); it++) {
		int size = *it;
		csd[size]++;
	}
}

template<typename number>
string VPCells<number>::get_coloured_mgl(number diameter) {
	string ret(Utils::sformat(".Box:%lf,%lf,%lf\n", this->_box->box_sides().x, this->_box->box_sides().y, this->_box->box_sides().z));

	number radius = diameter/2.;

	vector<LR_vector<number> > colors;
	colors.push_back(LR_vector<number>(0.f, 0.f, 1.f));
	colors.push_back(LR_vector<number>(0.f, 1.f, 0.f));
	colors.push_back(LR_vector<number>(1.f, 0.f, 0.f));
	colors.push_back(LR_vector<number>(1.f, 0.f, 1.f));
	colors.push_back(LR_vector<number>(1.f, 0.6f, 0.f));
	colors.push_back(LR_vector<number>(1.f, 1.f, 0.f));
	colors.push_back(LR_vector<number>(0.5f, 0.5f, 0.5f));
	colors.push_back(LR_vector<number>(0.56f, 0.f, 1.f));
	colors.push_back(LR_vector<number>(0.647f, 0.165f, 0.165f));
	colors.push_back(LR_vector<number>(0.f, 1.f, 1.f));
	colors.push_back(LR_vector<number>(0.98f, 0.855f, 0.867f));
	typename vector<LR_vector<number> >::iterator curr_color_it = colors.begin();

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		int c_id = _clusters[p->index];
		if(_color_map.find(c_id) == _color_map.end()) {
			_color_map[c_id] = *curr_color_it;
			curr_color_it++;
			if(curr_color_it == colors.end()) curr_color_it = colors.begin();
		}

		LR_vector<number> color = _color_map[c_id];
		ret += Utils::sformat("%lf %lf %lf @ %lf C[%lf,%lf,%lf]\n", p->pos.x, p->pos.y, p->pos.z, radius, color.x, color.y, color.z);
	}

	return ret;
}

template class VPCells<float>;
template class VPCells<double>;

template<typename number>
VoidPercolation<number>::VoidPercolation() {
	_probe_diameter = _particle_diameter = 1.;
	_cells = NULL;
	_probe_cells = _replica_probe_cells = NULL;
	_probes = _replica_probes = NULL;
	_probe = NULL;
	_particles = NULL;
	_insertions = _N = _N_replica = 0;
	_replica_box = NULL;
	_print_mgl = false;
}

template<typename number>
VoidPercolation<number>::~VoidPercolation() {
	if(_cells != NULL) delete _cells;
	if(_probe_cells != NULL) delete _probe_cells;
	if(_replica_probe_cells != NULL) delete _replica_probe_cells;
	if(_probes != NULL) {
		for(int i = 0; i < _insertions; i++) delete _probes[i];
		delete _probes;
	}
	if(_replica_probes != NULL) {
		for(int i = 0; i < 8*_insertions; i++) delete _replica_probes[i];
		delete _replica_probes;
	}
	if(_probe != NULL) delete _probe;
	if(_particles != NULL) delete _particles;
	if(_replica_box != NULL) delete _replica_box;
}

template<typename number>
void VoidPercolation<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "insertions", &_insertions, 1);
	getInputNumber(&my_inp, "probe_diameter", &_probe_diameter, 1);
	getInputNumber(&my_inp, "particle_diameter", &_particle_diameter, 0);
	getInputBool(&my_inp, "print_mgl", &_print_mgl, 0);

	_replica_box = BoxFactory::make_box<number>(sim_inp);
	_replica_box->get_settings(sim_inp);
}

template<typename number>
void VoidPercolation<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	_N = *config_info.N + 1;
	_N_replica = 8*_insertions;

	_probes = new BaseParticle<number>*[_insertions];
	_replica_probes = new BaseParticle<number>*[_N_replica];
	for(int i = 0; i < _insertions; i++) {
		_probes[i] = new BaseParticle<number>();
		BaseParticle<number> *p = _probes[i];
		p->init();
		p->index = i;
		p->type = P_A;
		p->pos = LR_vector<number>(0., 0., 0.);

		for(int j = 0; j < 8; j++) {
			int index = i*8 + j;
			_replica_probes[index] = new BaseParticle<number>();
			BaseParticle<number> *q = _replica_probes[index];
			q->init();
			q->index = index;
			q->type = P_A;
			q->pos = LR_vector<number>(0., 0., 0.);
		}
	}

	_probe = new BaseParticle<number>();
	_probe->init();
	_probe->index = _N - 1;
	_probe->type = P_A;

	// we make a copy of the _particles array and add the probe as an additional particle at the end of it
	_particles = new BaseParticle<number> *[_N];
	for(int i = 0; i < _N-1; i++) _particles[i] = config_info.particles[i];
	_particles[_N - 1] = _probe;

	_rcut = (_probe_diameter + _particle_diameter)/2.;
	_sqr_rcut = SQR(_rcut);
	_cells = new Cells<number>(_N, config_info.box);
	_cells->init(_particles, _rcut);

	_probe_cells = new VPCells<number>(_insertions, config_info.box);
	_probe_cells->init(_probes, _probe_diameter);

	LR_vector<number> box = config_info.box->box_sides();
	_replica_box->init(2.*box.x, 2.*box.y, 2.*box.z);
	_replica_probe_cells = new VPCells<number>(_N_replica, _replica_box);
	_replica_probe_cells->init(_replica_probes, _probe_diameter);

	if(_insertions == 0) {
		number tot_V = config_info.box->box_sides().x*config_info.box->box_sides().y*config_info.box->box_sides().z;
		number probe_V = M_PI*pow(_probe_diameter, 3.)/6.;
		_insertions = (int)(tot_V/probe_V)*50;
		OX_LOG(Logger::LOG_INFO, "VoidPercolation: insertions = %d", _insertions);
	}
}

template<typename number>
string VoidPercolation<number>::get_output_string(llint curr_step) {
	stringstream outstr;
	LR_vector<number> box = this->_config_info.box->box_sides();

	_cells->global_update();
	for(int i = 0; i < _insertions; i++) {
		bool overlap;
		int tries = 0;
		do {
			overlap = false;
			_probe->pos = LR_vector<number> (drand48()*box.x, drand48()*box.y, drand48()*box.z);
			_cells->single_update(_probe);

			vector<BaseParticle<number> *> neighs = _cells->get_complete_neigh_list(_probe);
			typename vector<BaseParticle<number> *>::iterator it;
			for(it = neighs.begin(); it != neighs.end() && !overlap; it++) {
				BaseParticle<number> *q = *it;
				LR_vector<number> dr = this->_config_info.box->min_image(_probe, q);
				if(dr.norm() < _sqr_rcut) overlap = true;
			}
			tries++;
			if(tries > 10000000) throw oxDNAException("Can't insert ghost particle n. %d", i);
		} while(overlap);

		_probes[i]->pos = _probe->pos;
		for(int x = 0; x < 2; x++) {
			LR_vector<number> shift(0., 0., 0.);
			shift[0] = x*box.x;
			for(int y = 0; y < 2; y++) {
				shift[1] = y*box.y;
				for(int z = 0; z < 2; z++) {
					shift[2] = z*box.z;

					int idx = 8*i + x + 2*(y + 2*z);
					_replica_probes[idx]->pos = _probe->pos + shift;
				}
			}
		}
	}

	_probe_cells->global_update(true);
	_probe_cells->compute_csd();

	if(_print_mgl) return _probe_cells->get_coloured_mgl(_probe_diameter);
	else {
		_replica_probe_cells->global_update(true);
		_replica_probe_cells->compute_csd();

//		for(int i = 1; i <= _insertions; i++) {
//			if(_probe_cells->csd[i] > 0) {
//				if(_probe_cells->csd[i] != 8*_replica_probe_cells->csd[i]) {
//					outstr << "PERCOLATION " <<
//					printf("PERCOLATION%d %d\n", i, _probe_cells->csd[i]);
//				}
//			}
//		}
		bool percolated = false;
		for(int i = _insertions + 1; i <= _N_replica; i++) {
			int n = _replica_probe_cells->csd[i];
			if(n > 0) {
				percolated = true;
				outstr << "1 PERCOLATION IN ";
				// check if we percolate in 1, 2 or 3 dimensions
				if(i/2 <= _insertions && _probe_cells->csd[i/2] == n/4) outstr << "ONE DIMENSION" << endl;
				else if(i/4 <= _insertions && _probe_cells->csd[i/4] == n/2) outstr << "TWO DIMENSIONS" << endl;
				else if(i/8 <= _insertions && _probe_cells->csd[i/8] == n) outstr << "THREE DIMENSIONS" << endl;
				else outstr << Utils::sformat("CAN'T FIND OUT, SOMETHING'S WRONG (%d %d)\n", i, n) << endl;
			}
		}
		if(!percolated) outstr << "0" << endl;
	}

//	outstr << n_overlaps / (number)_insertions << " " << n_overlaps;

	return outstr.str();
}

template class VoidPercolation<float>;
template class VoidPercolation<double>;

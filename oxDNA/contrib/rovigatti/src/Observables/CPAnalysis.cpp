/*
 * CPAnalysis.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "CPAnalysis.h"
#include <sstream>

template<typename number>
CPCells<number>::CPCells(int &N, BaseBox<number> *box) : Cells<number>(N, box) {

}

template<typename number>
CPCells<number>::~CPCells() {

}

template<typename number>
number CPCells<number>::assign_density_to_cells() {
	if(this->_N_cells_side[0] != this->_N_cells_side[1] || this->_N_cells_side[0] != this->_N_cells_side[2]) throw oxDNAException("CPAnalysis does not support non-cubic boxes");
	if(this->_N_cells_side[0] < 5) throw oxDNAException("The system should be large enough to contain 5x5x5 cells");

	number current_sigma = this->_box_sides[0]/this->_N_cells_side[0];
	number max_dist = 2*current_sigma;
	number sqr_max_dist = SQR(max_dist);

	cell_density.resize(this->_N_cells, 0);
	coarse_grained_cell_density.resize(this->_N_cells, 0);

	// we start by assign a density for each cell
	for(int c_ind = 0; c_ind < this->_N_cells; c_ind++) {
		int ind[3] = {
			c_ind % this->_N_cells_side[0],
			(c_ind / this->_N_cells_side[0]) % this->_N_cells_side[1],
			c_ind / (this->_N_cells_side[0]*this->_N_cells_side[1])
		};

		LR_vector<number> cell_pos((ind[0] + 0.5)*current_sigma, (ind[1] + 0.5)*current_sigma, (ind[2] + 0.5)*current_sigma);
		int N = 0;

		int loop_ind[3];
		for(int j = -2; j < 3; j++) {
			loop_ind[0] = (ind[0] + j + this->_N_cells_side[0]) % this->_N_cells_side[0];
			for(int k = -2; k < 3; k++) {
				loop_ind[1] = (ind[1] + k + this->_N_cells_side[1]) % this->_N_cells_side[1];
				for(int l = -2; l < 3; l++) {
					loop_ind[2] = (ind[2] + l + this->_N_cells_side[2]) % this->_N_cells_side[2];
					int loop_index = loop_ind[0] + this->_N_cells_side[0]*(loop_ind[1] + loop_ind[2]*this->_N_cells_side[1]);

					BaseParticle<number> *q = this->_heads[loop_index];
					while(q != P_VIRTUAL) {
						if(this->_box->sqr_min_image_distance(q->pos, cell_pos) < sqr_max_dist) N++;
						q = this->_next[q->index];
					}
				}
			}
		}

		cell_density[c_ind] = N/(4.*M_PI*CUB(max_dist)/3.);
	}

	// and then we smooth it over its first neighbours
	for(int c_ind = 0; c_ind < this->_N_cells; c_ind++) {
		int ind[3] = {
			c_ind % this->_N_cells_side[0],
			(c_ind / this->_N_cells_side[0]) % this->_N_cells_side[1],
			c_ind / (this->_N_cells_side[0]*this->_N_cells_side[1])
		};

		coarse_grained_cell_density[c_ind] = 2*cell_density[c_ind];

		// for each direction
		for(int i = 0; i < 3; i++) {
			for(int shift = -1; shift < 2; shift += 2) {
				int neigh_ind[3] = {ind[0], ind[1], ind[2]};
				neigh_ind[i] = (neigh_ind[i] + shift + this->_N_cells_side[i]) % this->_N_cells_side[i];
				int neigh_index = neigh_ind[0] + this->_N_cells_side[0]*(neigh_ind[1] + neigh_ind[2]*this->_N_cells_side[1]);
				coarse_grained_cell_density[c_ind] += cell_density[neigh_index];
			}
		}
	}

	// and we normalise all the densities
	for(int c_ind = 0; c_ind < this->_N_cells; c_ind++) {
		coarse_grained_cell_density[c_ind] /= 8.;
	}

	return current_sigma;
}

template<typename number>
int CPCells<number>::get_N_particles_in_cell(int cell_index) {
	int N = 0;
	BaseParticle<number> *q = this->_heads[cell_index];
	while(q != P_VIRTUAL) {
		q = this->_next[q->index];
		N++;
	}

	return N;
}

template class CPCells<float>;
template class CPCells<double>;

template<typename number>
CPAnalysis<number>::CPAnalysis(): BaseObservable<number>() {
	_type = -1;
	_cells = NULL;
}

template<typename number>
CPAnalysis<number>::~CPAnalysis() {
	if(_cells != NULL) delete _cells;
}

template<typename number>
void CPAnalysis<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "particle_type", &_type, 0);
	getInputNumber(&my_inp, "sigma", &_sigma, 1);
}

template<typename number>
void CPAnalysis<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	_cells = new CPCells<number>(*config_info.N, config_info.box);
	if(_type != -1) _cells->set_allowed_type(_type);
	_cells->init(config_info.particles, _sigma);
}

template<typename number>
std::string CPAnalysis<number>::get_output_string(llint curr_step) {
	_cells->global_update(true);
	number current_sigma = _cells->assign_density_to_cells();

	stringstream ss;
	for(int i = 0; i < _cells->get_N_cells(); i++) ss << _cells->coarse_grained_cell_density[i] << " " << _cells->cell_density[i] << " " << current_sigma << endl;

	return ss.str();
}

template class CPAnalysis<float>;
template class CPAnalysis<double>;

/*
 * SRDThermostat.cpp
 *
 *  Created on: 05/dic/2013
 *      Author: lorenzo
 */

#include <cfloat>

#include "SRDThermostat.h"

template<typename number>
SRDThermostat<number>::SRDThermostat(number &box_side) : BaseThermostat<number>(), _box_side(box_side) {
	_cells = NULL;
	_srd_particles = NULL;

	_N_particles = 0;
	_N_cells = _N_cells_side = _N_per_cell = 0;
	_apply_every = 0;
}

template<typename number>
SRDThermostat<number>::~SRDThermostat() {
	if(_cells != NULL) delete[] _cells;
	if(_srd_particles != NULL) delete[] _srd_particles;
}

template<typename number>
void SRDThermostat<number>::get_settings(input_file &inp) {
	BaseThermostat<number>::get_settings(inp);

	float tmp;
	getInputFloat(&inp, "srd_r_cell", &tmp, 1);
	_r_cell = (number) tmp;

	getInputFloat(&inp, "dt", &tmp, 1);
	_dt = (number) tmp;

	getInputFloat(&inp, "srd_alpha", &tmp, 1);
	_alpha = (number) tmp;

	getInputInt(&inp, "srd_N_per_cell", &_N_per_cell, 1);
	getInputInt(&inp, "srd_apply_every", &_apply_every, 1);

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number>(raw_T);
}

template<typename number>
void SRDThermostat<number>::init(int N_part) {
	BaseThermostat<number>::init(N_part);

	_N_cells_side = (int) (floor(_box_side / _r_cell) + 0.1); 
	_N_cells = _N_cells_side * _N_cells_side * _N_cells_side;
	_r_cell = _box_side / _N_cells_side;
	_cells = new SRDCell<number>[_N_cells];

	_N_particles = _N_cells * _N_per_cell;
	_srd_particles = new SRDParticle<number>[_N_particles + N_part];

	// the sum of the srd particle masses in a cell should be equal
	// to the mass of a particle (which is fixed to 1)
	_m = 1. / _N_per_cell;

	number rescale_factor = sqrt(_T/_m);
	for(int i = 0; i < _N_particles; i++) {
		_srd_particles[i].r.x = drand48() * _box_side;
		_srd_particles[i].r.y = drand48() * _box_side;
		_srd_particles[i].r.z = drand48() * _box_side;

		_srd_particles[i].v.x = Utils::gaussian<number>() * rescale_factor;
		_srd_particles[i].v.y = Utils::gaussian<number>() * rescale_factor;
		_srd_particles[i].v.z = Utils::gaussian<number>() * rescale_factor;
	}
}

template<typename number>
int SRDThermostat<number>::_get_cell_index(LR_vector<number> &r) {
	int ind[3];
	ind[0] = (int) ((r.x / _box_side - floor(r.x / _box_side)) * (1. - DBL_EPSILON) * _N_cells_side);
	ind[1] = (int) ((r.y / _box_side - floor(r.y / _box_side)) * (1. - DBL_EPSILON) * _N_cells_side);
	ind[2] = (int) ((r.z / _box_side - floor(r.z / _box_side)) * (1. - DBL_EPSILON) * _N_cells_side);

	return (ind[0] * _N_cells_side + ind[1]) * _N_cells_side + ind[2];
}

template<typename number>
void SRDThermostat<number>::apply(BaseParticle<number> **particles, llint curr_step) {
	if(curr_step % _apply_every == 0) {
		for(int i = 0; i < _N_cells; i++) {
			_cells[i].v_avg = LR_vector<number>(0., 0., 0.);
			_cells[i].tot_mass = 0.;
			_cells[i].head = NULL;
		}

		number kT = 0.f;
		for(int i = 0; i < _N_particles; i++) {
			SRDParticle<number> *p = _srd_particles + i;
			p->r += p->v * _dt * _apply_every;

			int ind = _get_cell_index(p->r);
			SRDCell<number> *c = _cells + ind;
			c->v_avg += _m * p->v;
			c->tot_mass += _m;

			p->next = c->head;
			c->head = p;

			kT += p->v.norm() * _m;
		}

		for(int i = 0; i < this->_N_part; i++) {
			BaseParticle<number> *p = particles[i];

			int ind = _get_cell_index(p->pos);
			SRDCell<number> *c = _cells + ind;
			c->v_avg += p->vel;
			// regular particles have mass 1
			c->tot_mass += 1.;
			// regular particles have to be given a cell, so we use srd particles
			// to assign cell indexes and copy the velocity;
			SRDParticle<number> *fake_p = _srd_particles + (_N_particles + i);
			fake_p->v = p->vel;
			fake_p->next = c->head;
			c->head = fake_p;
		}

		for(int i = 0; i < _N_cells; i++) {
			SRDCell<number> *c = _cells + i;
			c->v_avg /= c->tot_mass;

			// rotation matrix
			LR_matrix<number> R = Utils::get_random_rotation_matrix_from_angle(_alpha);
			SRDParticle<number> *p = c->head;
			while(p != NULL) {
				p->v = R*(p->v - c->v_avg) + c->v_avg;
				p = p->next;
			}
		}

		// now we have to copy back the updated velocities of regular particles
		for(int i = 0; i < this->_N_part; i++) {
			BaseParticle<number> *p = particles[i];
			SRDParticle<number> *fake_p = _srd_particles + (_N_particles + i);
			p->vel = fake_p->v;
		}
	}
}

template class SRDThermostat<float>;
template class SRDThermostat<double>;

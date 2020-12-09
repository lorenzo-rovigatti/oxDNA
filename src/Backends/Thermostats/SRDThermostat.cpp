/*
 * SRDThermostat.cpp
 *
 *  Created on: 05/dic/2013
 *      Author: lorenzo/Flavio
 */

#include "SRDThermostat.h"

#include "../../Boxes/BaseBox.h"
#include "../../Utilities/ConfigInfo.h"

#include <cfloat>

SRDThermostat::SRDThermostat(BaseBox *box) :
				BaseThermostat(),
				_box(box) {
	_cells = nullptr;
	_srd_particles = nullptr;
	_N_particles = 0;
	_N_cells = _N_cells_side = _N_per_cell = 0;
	_apply_every = 0;
	_rescale_factor = 0.;
	_m = -1.;
	_r_cell = -1.;
	_is_cuda = false;
	_dt = -1.;
}

SRDThermostat::~SRDThermostat() {
	if(!_is_cuda) {
		if(_cells != nullptr) {
			delete[] _cells;
		}
		if(_srd_particles != nullptr) {
			delete[] _srd_particles;
		}
	}
}

void SRDThermostat::get_settings(input_file &inp) {
	BaseThermostat::get_settings(inp);

	CHECK_BOX("SRDThermostat", inp);

	float tmp;
	getInputFloat(&inp, "srd_r_cell", &tmp, 1);
	_r_cell = (number) tmp;

	getInputFloat(&inp, "dt", &tmp, 1);
	_dt = (number) tmp;

	getInputInt(&inp, "srd_N_per_cell", &_N_per_cell, 1);
	getInputInt(&inp, "srd_apply_every", &_apply_every, 1);

	char backend[512];
	getInputString(&inp, "backend", backend, 1);
	_is_cuda = (strcmp(backend, "CUDA") == 0);

	OX_LOG(Logger::LOG_INFO, "SRD thermostat: T=%g, r_cell: %g, N_per_cell: %d, rescale_factor: %g", _T, _r_cell, _N_per_cell, _rescale_factor);
}

void SRDThermostat::init() {
	BaseThermostat::init();

	_rescale_factor = sqrt(_T);

	number L = _box->box_sides()[0];
	_N_cells_side = (int) (floor(L / _r_cell) + 0.1);
	_N_cells = _N_cells_side * _N_cells_side * _N_cells_side;
	_r_cell = L / _N_cells_side;

	_N_particles = _N_cells * _N_per_cell;
	// the sum of the srd particle masses in a cell should be equal
	// to the mass of a particle (which is fixed to 1)
	_m = 1. / _N_per_cell;

	// we allocate memory only if we simulate on the CPU
	if(!_is_cuda && _cells == nullptr) {
		_cells = new SRDCell[_N_cells];
		_srd_particles = new SRDParticle[_N_particles + CONFIG_INFO->N()];

		number rescale_factor = sqrt(_T / _m);
		for(int i = 0; i < _N_particles; i++) {
			_srd_particles[i].r.x = drand48() * L;
			_srd_particles[i].r.y = drand48() * L;
			_srd_particles[i].r.z = drand48() * L;

			_srd_particles[i].v.x = Utils::gaussian() * rescale_factor;
			_srd_particles[i].v.y = Utils::gaussian() * rescale_factor;
			_srd_particles[i].v.z = Utils::gaussian() * rescale_factor;

			_srd_particles[i].L = LR_vector(0., 0., 0.);
		}
	}
}

int SRDThermostat::_get_cell_index(LR_vector &r) {
	LR_vector norm_in_box = _box->normalised_in_box(r) * _N_cells_side;
	int ind[3];
	ind[0] = (int) norm_in_box.x;
	ind[1] = (int) norm_in_box.y;
	ind[2] = (int) norm_in_box.z;

	return (ind[0] * _N_cells_side + ind[1]) * _N_cells_side + ind[2];
}

void SRDThermostat::apply(std::vector<BaseParticle *> &particles, llint curr_step) {
	if(_is_cuda) throw oxDNAException("The apply method of the SRD thermostat has been called on the CPU on a CUDA-enabled simulation. This should not happen.");

	apply1(particles, curr_step);
	//return apply2 (particles, curr_step);
	//return apply3 (particles, curr_step);
}

// Andersen-MPCD
// conserves linear momentum, does not conserve angular momentum
// angular momenta of the solute particles are refreshed
void SRDThermostat::apply1(std::vector<BaseParticle *> &particles, llint curr_step) {
	if(!(curr_step % _apply_every == 0)) return;

	number L = _box->box_sides()[0];
	for(int i = 0; i < _N_cells; i++) {
		_cells[i].P = LR_vector(0., 0., 0.);  // total momentum
		_cells[i].PR = LR_vector(0., 0., 0.); // sum of random components
		_cells[i].tot_mass = (number) 0.f;
		_cells[i].head = nullptr;
		int x, y, z;
		x = i % _N_cells_side;
		y = (i / _N_cells_side) % _N_cells_side;
		z = (i / _N_cells_side / _N_cells_side) % _N_cells_side;
		// center of cell
		_cells[i].centre = LR_vector((x + 0.5) * L / _N_cells_side, (y + 0.5) * L / _N_cells_side, (z + 0.5) * L / _N_cells_side);
	}

	// put particles in cells and find out average
	// velocity and angular momentum 
	number sqrt_m = sqrt(_m);
	for(int i = 0; i < _N_particles; i++) {
		SRDParticle *p = _srd_particles + i;
		p->r += p->v * _dt * _apply_every;

		int ind = _get_cell_index(p->r);
		p->cell_index = ind;
		SRDCell *c = _cells + ind;

		// need to add this BEFORE refreshing
		c->P += _m * p->v; // linear momentum 

		// new velocity; later on we correct to get conservation
		// of linar momentum within each cell. For now, we just
		// refresh the velocity completely
		p->v = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * (_rescale_factor / sqrt_m);

		// store the contribution to mass, linear and angular momentum of cell
		c->PR += _m * p->v;
		c->tot_mass += _m;

		p->next = c->head;
		c->head = p;
	}

	for(uint i = 0; i < particles.size(); i++) {
		BaseParticle *p = particles[i];
		int ind = _get_cell_index(p->pos);
		SRDCell *c = _cells + ind;

		// need to add this BEFORE refreshing
		c->P += p->vel; // linear momentum, mass = 1 

		// we refresh the linear and angular momentum completely; later on, we
		// will correct to get conservation of linear and angular momentum
		// within each cell
		p->vel = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
		p->L = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;

		// store the contribution to the total mass and spin inertia
		c->PR += p->vel;
		c->tot_mass += (number) 1.f;

		// regular particles have to be given a cell, so we use srd particles
		// to assign cell indexes and copy the velocity;
		SRDParticle *fake_p = _srd_particles + (_N_particles + i);
		fake_p->v = p->vel;
		fake_p->L = p->L;
		fake_p->cell_index = ind;

		fake_p->next = c->head;
		c->head = fake_p;
	}

	for(int i = 0; i < _N_cells; i++) {
		SRDCell *c = _cells + i;

		if(c->tot_mass > 0.) {
			c->P /= c->tot_mass;
			c->PR /= c->tot_mass;
		}

		//printf ("%d (%g %g) %g %g %g -- %g %g %g %g@@\n", i, c->tot_mass, c->tot_I, c->P.module(), c->L.module(), c->L_spin.module(), c->PR.module(), c->LR.module(), c->LR_spin.module(), c->dLgn.module()); 
	}

	// update velocities of SRD particles
	for(int i = 0; i < _N_particles; i++) {
		SRDParticle *p = _srd_particles + i;
		p->v += _cells[p->cell_index].P - _cells[p->cell_index].PR;
	}

	// update velocities and angular momenta
	// of solute particles
	//return;
	for(uint i = 0; i < particles.size(); i++) {
		BaseParticle *p = particles[i];
		SRDParticle *fake_p = _srd_particles + (_N_particles + i);
		SRDCell *c = _cells + fake_p->cell_index;
		p->vel = fake_p->v + c->P - c->PR;
		p->L = fake_p->L;
	}
}

// non working at the moment... 
//
//void SRDThermostat::apply2(std::vector<BaseParticle *> &particles, llint curr_step) {
//	if (!(curr_step % _apply_every == 0)) return;
//
//	for(int i = 0; i < _N_cells; i++) {
//		_cells[i].P = LR_vector(0., 0., 0.);  // total momentum
//		_cells[i].L = LR_vector(0., 0., 0.);  // total anguler momentum
//		_cells[i].L_spin = LR_vector(0., 0., 0.);  // total anguler momentum
//		_cells[i].PR = LR_vector(0., 0., 0.); // sum of random components
//		_cells[i].LR = LR_vector(0., 0., 0.); // sum of random components
//		_cells[i].LR_spin = LR_vector(0., 0., 0.); // sum of random components
//		_cells[i].dLgn = LR_vector (0., 0., 0.); // spin-orbit coupling
//		_cells[i].tot_mass = (number) 0.f;
//		_cells[i].tot_I = (number) 0.f;
//		_cells[i].head = NULL;
//		int x, y, z;
//		x = i % _N_cells_side;
//		y = (i / _N_cells_side) % _N_cells_side;
//		z = (i / _N_cells_side / _N_cells_side) % _N_cells_side;
//		_cells[i].centre = LR_vector ((x + 0.5) * _box_side / _N_cells_side,
//											  (y + 0.5) * _box_side / _N_cells_side,
//											  (z + 0.5) * _box_side / _N_cells_side); // center of cell
//	}
//
//	// put particles in cells and find out average
//	// velocity and angular momentum
//	number sqrt_m = sqrt(_m);
//	for(int i = 0; i < _N_particles; i++) {
//		SRDParticle *p = _srd_particles + i;
//		p->r += p->v * _dt * _apply_every;
//
//		int ind = _get_cell_index(p->r);
//		p->cell_index = ind;
//		SRDCell *c = _cells + ind;
//
//		// need to add this BEFORE refreshing
//		c->P += _m * p->v; // linear momentum
//		c->L += _m * p->r.minimum_image (c->centre, _box_side).cross(p->v); // orbit term
//
//		// new velocity; later on we correct to get conservation
//		// of linar momentum within each cell. For now, we just
//		// refresh the velocity completely
//		p->v = LR_vector (Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor / sqrt_m;
//
//		// store the contribution to mass, linear and angular momentum of cell
//		c->PR += _m * p->v;
//		c->LR += _m * p->r.minimum_image(c->centre, _box_side).cross(p->v);
//		c->tot_mass += _m;
//
//		p->next = c->head;
//		c->head = p;
//	}
//
//	for(int i = 0; i < this->_N_part; i++) {
//		BaseParticle *p = particles[i];
//
//		int ind = _get_cell_index(p->pos);
//		SRDCell *c = _cells + ind;
//
//		// need to add this BEFORE refreshing
//		c->P += p->vel; // linear momentum, mass = 1
//		c->L += p->pos.minimum_image (c->centre, _box_side).cross (p->vel); // orbit term, mass = 1
//		c->L_spin += p->L; // spin term, inertia = 1
//
//		// we refresh the linear and angular momentum completely; later on, we
//		// will correct to get conservation of linear and angular momentum
//		// within each cell
//		p->vel = LR_vector (Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
//		p->L = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _rescale_factor;
//
//		// store the contribution to the total mass and spin inertia
//		c->PR += p->vel;
//		c->LR += p->pos.minimum_image (c->centre, _box_side).cross (p->vel); // orbit term, mass = 1
//		c->LR_spin += p->L;
//		c->tot_mass += (number) 1.f;
//		c->tot_I += (number) 1.f;
//
//		// regular particles have to be given a cell, so we use srd particles
//		// to assign cell indexes and copy the velocity;
//		SRDParticle *fake_p = _srd_particles + (_N_particles + i);
//		fake_p->v = p->vel;
//		fake_p->L = p->L;
//		fake_p->cell_index = ind;
//
//		fake_p->next = c->head;
//		c->head = fake_p;
//	}
//
//	for(int i = 0; i < _N_cells; i++) {
//		SRDCell *c = _cells + i;
//
//		if (c-> tot_mass > 0.) {
//			c->P /= c->tot_mass;
//			c->PR /= c->tot_mass;
//		}
//
//		if (c->tot_I > 0.) {
//			c->dLgn = (c->LR - c->L) / c->tot_I; // new minus old angular momentum, to be balanced later
//			c->L /= c->tot_I;
//			c->L_spin /= c->tot_I;
//			c->LR /= c->tot_I;
//			c->LR_spin /= c->tot_I;
//		}
//		if (c->tot_I > 0.)
//		printf ("%d (%g %g) %g %g %g -- %g %g %g %g@@\n", i, c->tot_mass, c->tot_I, c->P.module(), c->L.module(), c->L_spin.module(), c->PR.module(), c->LR.module(), c->LR_spin.module(), c->dLgn.module());
//	}
//
//	// update velocities of SRD particles
//	for(int i = 0; i < _N_particles; i++) {
//		SRDParticle *p = _srd_particles + i;
//		p->v += _cells[p->cell_index].P - _cells[p->cell_index].PR;
//	}
//
//	// update velocities and angular momenta
//	// of solute particles
//	return;
//	for(int i = 0; i < this->_N_part; i++) {
//		BaseParticle *p = particles[i];
//		SRDParticle *fake_p = _srd_particles + (_N_particles + i);
//		SRDCell *c = _cells + fake_p->cell_index;
//		p->vel = fake_p->v + c->P - c->PR;
//		p->L = fake_p->L + c->L_spin - c->LR_spin - c->dLgn;
//	}
//}

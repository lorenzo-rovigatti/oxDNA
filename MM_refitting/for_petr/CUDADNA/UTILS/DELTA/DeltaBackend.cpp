/*
 * DeltaBackend.cpp
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

#include "DeltaBackend.h"
#include "../../IOManager.h"

#define DELTA 0.119

DeltaBackend::DeltaBackend(IOManager *IO) : MC_CPUBackend<double>(IO) {

}

DeltaBackend::~DeltaBackend() {
	for(int i = 0; i < 2; i++) delete _strands[i];
}

void DeltaBackend::init(ifstream &conf_input) {
	MCBackend<double>::init(conf_input);

	for(int i = 0; i < this->_N; i++) {
		// this is needed for the first _compute_energy()
		this->_particles[i].set_positions();
		this->_particles[i].orientationT = this->_particles[i].orientation.get_transpose();
	}

	if(_N != 48) _IO->die("Configuration must contain 48 nucleotides");

	// we have to bring back in the same box the two strands that form each cylinder
	Strand *mystrands[4];
	for(int i = 0; i < 4; i++) mystrands[i] = new Strand(_particles + i*12, 12);

	for(int i = 0; i < 4; i += 2) {
		LR_vector<double> disp = mystrands[i]->cm_pos - mystrands[i+1]->cm_pos;
		disp = LR_vector<double>(
					rint(disp.x / _box_side) * _box_side,
					rint(disp.y / _box_side) * _box_side,
					rint(disp.z / _box_side) * _box_side
				);
		mystrands[i+1]->translate(disp);
	}
	for(int i = 0; i < 4; i++) delete mystrands[i];

	// and now we create two strands, one for each cylinder
	for(int i = 0; i < 2; i++) _strands[i] = new Strand(_particles + i*24, 24);

	this->_update_lists();
	_compute_energy(false);
	_initial_U = _U;

	/*_patches_on0[0] = _strands[0]->first->pos + _strands[0]->first->pos_stack;
	_patches_on0[1] = _strands[0]->last->pos + _strands[0]->last->pos_stack;
	_patches_on0[2] = _strands[0]->first[11].pos + _strands[0]->first[11].pos_stack;
	_patches_on0[3] = _strands[0]->first[12].pos + _strands[0]->first[12].pos_stack;

	double ds[4] = {
		_strands[1]->cm_pos.sqr_distance(_strands[1]->first->pos + _strands[1]->first->pos_stack),
		_strands[1]->cm_pos.sqr_distance(_strands[1]->last->pos + _strands[1]->last->pos_stack),
		_strands[1]->cm_pos.sqr_distance(_strands[1]->first[11].pos + _strands[1]->first[11].pos_stack),
		_strands[1]->cm_pos.sqr_distance(_strands[1]->first[12].pos + _strands[1]->first[12].pos_stack)
	};

	_min_corona = 10000;
	_max_corona = 0;
	for(int i = 0; i < 4; i++) {
		if(ds[i] < _min_corona) _min_corona = ds[i];
		if(ds[1] > _max_corona) _max_corona = ds[i];
	}

	_min_corona = sqrt(_min_corona);
	_max_corona = sqrt(_max_corona) + CXST_RCHIGH;*/
}

void DeltaBackend::_compute_energy(bool all_interactions) {
	_U = (double) 0;
	for(int i = 0; i < _N; i++) {
		Particle<double> *p = &_particles[i];
		p->prepare_list();
		int neigh = p->next_neighbour();
		while(neigh != P_VIRTUAL) {
			// if all_interactions == false i and neigh interact if and only if they belong to the same cylinder
			if(neigh > i && (all_interactions || ((i / 24) == (neigh / 24)))) _U += _particle_particle_interaction(p, &_particles[neigh]);
			neigh = p->next_neighbour();
		}
	}
}

// see http://mathworld.wolfram.com/Sphere-SphereIntersection.html
double DeltaBackend::get_trial_volume() {
	double R = CXST_RCHIGH;
	double V = 4 * 4 * M_PI * R*R*R / 3.;

	LR_vector<double> dv = _strands[0]->first->pos +_strands[0]->first->pos_stack - (_strands[0]->last->pos +_strands[0]->last->pos_stack);
	double d = dv.module();
	if(d < 2*R) V -= M_PI * (4*R + d) * SQR(2*R - d) / 12.;

	dv = _strands[0]->first[11].pos + _strands[0]->first[11].pos_stack - (_strands[0]->first[12].pos + _strands[0]->first[12].pos_stack);
	d = dv.module();
	if(d < 2*R) V -= M_PI * (4*R + d) * SQR(2*R - d) / 12.;

	return V;
	//return 4 * M_PI * (pow(_max_corona, 3) - pow(_min_corona, 3)) / 3.;
}

/*double DeltaBackend::insertion_energy() {
	LR_vector<double> p11 = LR_vector<double>(0.5, 0, 0);
	LR_vector<double> p12 = LR_vector<double>(-0.5, 0, 0);

	LR_vector<double> new_pos = LR_vector<double>(2 * (1. + DELTA) * (drand48() - 0.5), 2 * (1. + DELTA) * (drand48() - 0.5), 2 * (1. + DELTA) * (drand48() - 0.5));
	// no overlap
	if(new_pos.norm() <= 1) return 0;

	LR_vector<double> patch = Utils::get_random_vector<double>();

	LR_vector<double> p21 = new_pos + 0.5 * patch;
	LR_vector<double> p22 = new_pos - 0.5 * patch;

	if(p21.distance(p11) <= DELTA) return -1;
	if(p22.distance(p11) <= DELTA) return -1;

	if(p21.distance(p12) <= DELTA) return -1;
	if(p22.distance(p12) <= DELTA) return -1;
	
	return 0;
	}*/

double DeltaBackend::insertion_energy() {
	int from[4] = {0, 11, 12, 23};
	// choose one of the 4 end-nucleotides on each cylinder
	Particle<double> *on0 = _strands[0]->first + from[lrand48() % 4];
	Particle<double> *on1 = _strands[1]->first + from[lrand48() % 4];

	LR_matrix<double> R = Utils::get_random_rotation_matrix<double>();
	_strands[1]->rotate(R);

	LR_vector<double> new_pos = Utils::get_random_vector_in_sphere<double>(CXST_RCHIGH);
	LR_vector<double> disp = on0->pos - on1->pos + new_pos + on0->pos_stack - on1->pos_stack;
	_strands[1]->translate(disp);

	this->_update_lists();
	this->_compute_energy(true);

	return _U - _initial_U;
}

/*double DeltaBackend::insertion_energy() {
	int from[4] = {0, 11, 12, 23};
	// choose one of the 4 end-nucleotides on each cylinder
	Particle<double> *on0 = _strands[0]->first + from[lrand48() % 4];

	LR_matrix<double> R = Utils::get_random_rotation_matrix<double>();
	_strands[1]->rotate(R);

	LR_vector<double> new_pos = Utils::get_random_vector_in_sphere<double>(_max_corona);
	while(new_pos.norm() < SQR(_min_corona)) new_pos = Utils::get_random_vector_in_sphere<double>(_max_corona);

	LR_vector<double> disp = on0->pos + on0->pos_stack - _strands[1]->cm_pos + new_pos;
	_strands[1]->translate(disp);

	// fast check on radial cut-off
	LR_vector<double> patches_on1[4] = {
		_strands[1]->first->pos + _strands[1]->first->pos_stack,
		_strands[1]->last->pos + _strands[1]->last->pos_stack,
		_strands[1]->first[11].pos + _strands[1]->first[11].pos_stack,
		_strands[1]->first[12].pos + _strands[1]->first[12].pos_stack
	};

  	double R2 = SQR(CXST_RCHIGH);
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			if(_patches_on0[i].sqr_distance(patches_on1[j]) < R2) {
				this->_update_lists();
				this->_compute_energy(true);
				
				return _U - _initial_U;
			}
		}
	}

	
	
	return 0;
}
*/

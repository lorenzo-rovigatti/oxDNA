/*
 * Particle.cpp
 *
 *  Created on: 21/set/2010
 *      Author: lorenzo
 */

#include "Particle.h"

template<typename number>
Particle<number>::Particle() : _N_neigh(0), _verlet_list(NULL), _principal_axis(1,0,0), _stack_axis(0,0,1), _third_axis(0,1,0), _N_ext_forces(0), index(-1), type(P_VIRTUAL) {
	en3 = (number) 0;
	en5 = (number) 0;
	esn3 = (number) 0;
	esn5 = (number) 0;
	inclust = false;
	e_neigh = 0;
	h_bonds = 0;
	ext_potential = (number) 0.;
	//_verlet_list = 0;
}

template<typename number>
void Particle<number>::copy_from(const Particle<number> &p) {
	index = p.index;
	type = p.type;
	btype = p.btype;
	_max_neigh = p._max_neigh;
	_N_neigh = p._N_neigh;
	pos = p.pos;
	vel = p.vel;
	orientation = p.orientation;
	orientationT = p.orientationT;
	pos_list = p.pos_list;
	force = p.force;
	en3 = p.en3;
	en5 = p.en5;
	esn3 = p.esn3;
	esn5 = p.esn5;
	n3 = p.n3;
	n5 = p.n5;
	pos_back = p.pos_back;
	pos_stack = p.pos_stack;
	pos_base = p.pos_base;
	ext_potential = p.ext_potential;

	//memcpy(_verlet_list, p._verlet_list, _N_neigh * sizeof(int));
	/*
	memcpy(_verlet_list, p._verlet_list, _max_neigh * sizeof(int));
	memcpy(e_neigh, p.e_neigh, (_max_neigh) * sizeof(number));
	memcpy(h_bonds, p.h_bonds, (_max_neigh) * sizeof(bool));
	*/
	memcpy(_verlet_list, p._verlet_list, _N_neigh * sizeof(int));
	memcpy(e_neigh, p.e_neigh, (_N_neigh) * sizeof(number));
	memcpy(h_bonds, p.h_bonds, (_N_neigh) * sizeof(bool));
}

template<typename number>
Particle<number>::~Particle() {
	if(_verlet_list != NULL) delete[] _verlet_list;
	if(e_neigh != NULL) delete[] e_neigh;
	if(h_bonds != NULL) delete[] h_bonds;
	for(int i = 0; i < _N_ext_forces; i++) delete _ext_forces[i];
}

template<typename number>
bool Particle<number>::add_ext_force(ExternalForce<number> *f) {
	if(_N_ext_forces == MAX_EXT_FORCES) return false;

	_ext_forces[_N_ext_forces] = f;
	_N_ext_forces++;

	return true;
}

template<typename number>
void Particle<number>::init(int max_neigh) {
	_max_neigh = max_neigh;
	_verlet_list = new int[_max_neigh];
	e_neigh = new number[_max_neigh];
	h_bonds = new bool[_max_neigh];

	_check();
}

template<typename number>
void Particle<number>::_check() {
	assert(index >= 0);
	//printf ("@@@ type %i\n", type);
	assert(type != P_VIRTUAL);
	assert(_verlet_list != NULL);
	assert(e_neigh != NULL);
	assert(h_bonds != NULL);
}

template<typename number>
void Particle<number>::_set_grooving(bool grooving) {
	_grooving = grooving;
}

template<typename number>
void Particle<number>::_set_back_position() {
	if (_grooving)
		pos_back = orientation*(_principal_axis*POS_BACK_A + _third_axis*POS_BACK_B);
	else{
		pos_back = orientation*_principal_axis*POS_BACK;
	}
}

template class Particle<double>;
template class Particle<float>;

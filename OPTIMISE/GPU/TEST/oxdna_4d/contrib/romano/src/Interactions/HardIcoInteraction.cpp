/*
 * HardIcoInteraction.cpp
 *
 *  Created on: 16/Lug/2018
 *      Author: Flavio 
 */

#include "HardIcoInteraction.h"

template<typename number>
HardIcoInteraction<number>::HardIcoInteraction() : BaseInteraction<number, HardIcoInteraction<number> >() {
	this->_int_map[HardIco] = &HardIcoInteraction<number>::_hi_pot;
	_close_vertexes = new int[12 * 5]; // 5 close vertexes for each vertex
	_tworinscribed = (number) (2.*sqrt((1./12.) + (1./(6.*sqrt(5.))))); // twice the radius of inscribed sphere
}

template<typename number>
HardIcoInteraction<number>::~HardIcoInteraction() {
	delete[] _close_vertexes;
}

template<typename number>
void HardIcoInteraction<number>::get_settings(input_file &inp) {
	BaseInteraction<number>::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512) && strncmp(tmps, "MC2", 512)) throw oxDNAException ("Cannot run Hard Icos with MD");
	
	float tmpf;
	getInputFloat (&inp, "patchy_delta", &tmpf, 1);
	_delta = (number) tmpf;
	
	OX_LOG(Logger::LOG_INFO, "Initializing HardIco interaction with delta %g", _delta);
	
	this->_rcut = (number) 1.001 * (1. + _delta);

	OX_LOG(Logger::LOG_INFO, "Using r_cut of %g", this->_rcut);
}

template<typename number>
void HardIcoInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);

	// compute the close vertexes thanks to a temporary hard icosahedron
	Icosahedron<number> tmp (6);
	tmp.pos = LR_vector<number> ((number)0.f, (number)0.f, (number)0.f);
	tmp.orientation = LR_matrix<number> ((number)1.f, (number)0.f, (number)0.f,
								 		(number)0.f, (number)1.f, (number)0.f,
								 		(number)0.f, (number)0.f, (number)1.f);
	tmp.orientationT = tmp.orientation.get_transpose();
	tmp.set_positions();

	for (int i = 0; i < 12; i ++) {
		int myn = 0;
		for (int j = 0; j < 12; j ++) {
			if (i != j and ((tmp.int_centers[i])*(tmp.int_centers[j]) > 0.)) {
				_close_vertexes[5*i + myn] = j;
				myn ++;
			}
		}
		if (myn != 5) throw oxDNAException("something wrong while initializing hardico interaction...");
	}

	// now we need to order the arrays so that we move counterclockwise
	for (int i = 0; i < 12; i ++) { // for each vertex
		bool good = false;
		while (good == false) {
			good = true;
			for (int j = 0; j < 4; j ++) {
				if (tmp.int_centers[_close_vertexes[5*i+j]]*tmp.int_centers[_close_vertexes[5*i + (j + 1)]] < 0.) {
					int save = _close_vertexes[5*i + j + 1];
					for (int k=1; k < 5-j-1; k ++) {
						_close_vertexes[5*i + j + k] = _close_vertexes[5*i + j + k + 1];
					}
					_close_vertexes[5*i+4] = save;
					good = false;
					break;
				}
			}
		}
	}

	// double check that all vectors are close to each other;
	bool good = true;
	for (int i = 0; i < 12; i ++) { // for each vertex
		for (int j = 0; j < 5; j ++) {
			number a = tmp.int_centers[i] * tmp.int_centers[_close_vertexes[5*i+j]];
			number b = tmp.int_centers[i] * tmp.int_centers[_close_vertexes[5*i+((j+1)%5)]];
			number c = tmp.int_centers[_close_vertexes[5*i+j]] * tmp.int_centers[_close_vertexes[5*i+((j+1)%5)]];
			if (a <= (number) 0.) good = false; 
			if (b <= (number) 0.) good = false; 
			if (c <= (number) 0.) good = false; 
		}
	}
	if (good == false) throw oxDNAException("not good at all, mfs");
}

template<typename number>
void HardIcoInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new Icosahedron<number>(6);
}

template<typename number>
number HardIcoInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

template<typename number>
number HardIcoInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number HardIcoInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _hi_pot (p, q, compute_r, update_forces);
}

template<typename number>
void HardIcoInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class HardIcoInteraction<float>;
template class HardIcoInteraction<double>;

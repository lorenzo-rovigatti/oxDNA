/*
 * BoxInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "BoxInteraction.h"
#include "InteractionUtils.h"

template<typename number>
BoxInteraction<number>::BoxInteraction() : BaseInteraction<number, BoxInteraction<number> >() {
	this->_int_map[Box] = &BoxInteraction<number>::_box_pot;
}

template<typename number>
BoxInteraction<number>::~BoxInteraction() {

}

template<typename number>
void BoxInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512)) throw oxDNAException ("Cannot run Box with MD");
	
	// let's get the three dimension
	getInputString (&inp, "box_sides", (char *)tmps, 1);
	float tmpf[3];
	
	int tmpi = sscanf (tmps, "%g, %g, %g", tmpf, tmpf + 1, tmpf + 2);
	_sides[0] = (number) tmpf[0];
	_sides[1] = (number) tmpf[1];
	_sides[2] = (number) tmpf[2];
	if (tmpi != 3) throw oxDNAException ("Could not read box sides from string \"%s\"", tmps);
	_lx = _sides[0];
	_ly = _sides[1];
	_lz = _sides[2];
	_largest = (_lx > _ly ? (_lx > _lz ? _lx : _lz) : (_ly > _lz ? _ly : _lz));
	_smallest = (_lx < _ly ? (_lx < _lz ? _lx : _lz) : (_ly < _lz ? _ly : _lz));
	
	OX_LOG(Logger::LOG_INFO, "Using box sides %g, %g, %g (%g, %g)", _sides[0], _sides[1], _sides[2], _smallest, _largest);
	
	this->_rcut = (number) 1.001 * sqrt(_sides[0] * _sides[0] + _sides[1] * _sides[1] + _sides[2] * _sides[2]);
	OX_LOG(Logger::LOG_INFO, "Using r_cut of %g", this->_rcut);
	//throw oxDNAException ("stopped here...");
}

template<typename number>
void BoxInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void BoxInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
number BoxInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number BoxInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number BoxInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _box_pot (p, q, r, update_forces);
}

template<typename number>
void BoxInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
bool BoxInteraction<number>::generate_random_configuration_overlap (BaseParticle<number> *p, BaseParticle<number> *q, number box_side) {
	LR_vector<number> dr = p->pos.minimum_image (q->pos, box_side);
	return InteractionUtils::box_overlap (p, q, dr, _lx, _ly, _lz);
}

template class BoxInteraction<float>;
template class BoxInteraction<double>;

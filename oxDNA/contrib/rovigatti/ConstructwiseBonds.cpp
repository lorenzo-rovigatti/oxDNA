/*
 * StrandwiseBonds.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Flavio Romano
 */

#include <sstream>

#include "ConstructwiseBonds.h"

using namespace std;

template<typename number>
ConstructwiseBonds<number>::ConstructwiseBonds() : _construct_strand_size(4), _construct_number(0) {

}

template<typename number>
ConstructwiseBonds<number>::~ConstructwiseBonds() {

}

template<typename number>
void ConstructwiseBonds<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "construct_size", &_construct_strand_size, 0);
}

template<typename number>
void ConstructwiseBonds<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	// we first count the number of constructs in the system so that we can
	// initialize the vector storing their centres of mass
	for(int i = 0; i < *config_info.N; i++) {
		BaseParticle<number> *p = config_info.particles[i];
		int p_construct_id = p->strand_id / _construct_strand_size;
		if(p_construct_id >= _construct_number) {
			_construct_number = p_construct_id+1;
			_construct_size.resize(_construct_number, 0);
			_construct_coms.resize(_construct_number, LR_vector<number>(0., 0., 0.));
		}
		_construct_size[p_construct_id]++;
	}
}

template<typename number>
string ConstructwiseBonds<number>::get_output_string(llint curr_step) {
	stringstream outstr;
	outstr << "# step " << curr_step << "\n";

	// compute the centres of mass of all the constructs
	for(int i = 0; i < _construct_number; i++) _construct_coms[i] = LR_vector<number>(0., 0., 0.);

	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		int p_construct_id = p->strand_id / _construct_strand_size;
		_construct_coms[p_construct_id] += p->get_abs_pos(*this->_config_info.box_side);
	}

	for(int i = 0; i < _construct_number; i++) _construct_coms[i] /= _construct_size[i];

	// compute the interaction map between constructs
	map<pair<int, int>, int> hbmap;
	vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();
	typename vector<ParticlePair<number> >::iterator it;
	for(it = pairs.begin(); it != pairs.end(); it ++ ) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;
		int p_construct_id = p->strand_id / _construct_strand_size;
		int q_construct_id = q->strand_id / _construct_strand_size;
		if(p_construct_id != q_construct_id) {
			number ene = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
			if (ene < HB_CUTOFF) {
				pair<int, int> k(p_construct_id, q_construct_id);
				if(hbmap.count(k) == 0) hbmap[k] = 1;
				else hbmap[k]++;
			}
		}
	}

	// print out everything
	map<pair<int, int>, int>::iterator it2;
	for(it2 = hbmap.begin(); it2 != hbmap.end(); it2 ++) {
		int i1, i2, nb;
		i1 = (*it2).first.first;
		i2 = (*it2).first.second;
		nb = (*it2).second;
		LR_vector<number> dist = _construct_coms[i1].minimum_image(_construct_coms[i2], *this->_config_info.box_side);

		outstr << i1 << " " << i2 << " " << nb << " " << dist.x << " " << dist.y << " " << dist.z << " " << dist.module() << "\n";

	}

	return outstr.str();
}

template class ConstructwiseBonds<float>;
template class ConstructwiseBonds<double>;

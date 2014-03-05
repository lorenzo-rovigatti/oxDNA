/*
 * StrandwiseBonds.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Flavio Romano
 */

#include "StrandwiseBonds.h"

#include <sstream>

template<typename number>
StrandwiseBonds<number>::StrandwiseBonds() {

}

template<typename number>
StrandwiseBonds<number>::~StrandwiseBonds() {

}

template<typename number>
std::string StrandwiseBonds<number>::get_output_string(llint curr_step) {
	std::stringstream outstr;
	outstr << "# step " << curr_step << "\n";

	std::map<std::pair<int, int>, int> hbmap;
	/*
	for (int i = 0; i < *this->_config_info.N; i++) {
		p = this->_config_info.particles[i];
		std::vector<BaseParticle<number> *> neighs = this->_config_info.interaction->get_neighbours(p, this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);
		typename std::vector<BaseParticle<number> *>::iterator it;
		for (it = neighs.begin(); it != neighs.end(); it ++) {
			BaseParticle<number> * q = (*it);
			if (p->strand_id < q->strand_id) {
				number ene = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
				if (ene < HB_CUTOFF) {
					std::pair<int, int> k (p->strand_id, q->strand_id);
					if (hbmap.count(k) == 0) hbmap[k] = 1;
					else hbmap[k] ++;
				}
			}
		}
	}
	*/
	std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> > pairs = this->_config_info.interaction->get_potential_interactions (this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);
	typename std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> >::iterator it;
	for (it = pairs.begin(); it != pairs.end(); it ++ ) {
		BaseParticle<number> * p = (*it).first;
		BaseParticle<number> * q = (*it).second;
		if (p->strand_id < q->strand_id) {
			number ene = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
			if (ene < HB_CUTOFF) {
				std::pair<int, int> k (p->strand_id, q->strand_id);
				if (hbmap.count(k) == 0) hbmap[k] = 1;
				else hbmap[k] ++;
			}
		}
	}

	std::map<std::pair<int, int>, int>::iterator it2;
	for (it2 = hbmap.begin(); it2 != hbmap.end(); it2 ++) {
		int i1, i2, nb;
		i1 = (*it2).first.first;
		i2 = (*it2).first.second;
		nb = (*it2).second;
		outstr << i1 << " " << i2 << " " << nb << "\n";
	}

	return outstr.str();
}

template class StrandwiseBonds<float>;
template class StrandwiseBonds<double>;

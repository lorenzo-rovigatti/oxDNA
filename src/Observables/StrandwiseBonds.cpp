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
	std::vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();
	typename std::vector<ParticlePair<number> >::iterator it;
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

/*
 * StressTensor.cpp
 *
 *  Created on: 11/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "StressTensor.h"

using namespace std;

template<typename number>
StressTensor<number>::StressTensor() : BaseObservable<number>() {
	_print_averaged_off_diagonal = false;
}

template<typename number>
StressTensor<number>::~StressTensor() {

}

template<typename number>
void StressTensor<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	int tmp = 0;
	getInputBoolAsInt(&my_inp, "avg_off_diagonal", &tmp, 0);
	_print_averaged_off_diagonal = (bool)tmp;
}

template<typename number>
string StressTensor<number>::get_output_string(llint curr_step) {
	vector<pair<BaseParticle<number> *, BaseParticle<number>*> > pairs = this->_config_info.interaction->get_potential_interactions(this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);

	vector<vector<number> > st(3, vector<number>(3, (number) 0));

	// now we loop on all the pairs in order to update the forces
	typename std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> >::iterator it;
	for (it = pairs.begin(); it != pairs.end(); it ++ ) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;
		p->force = LR_vector<number>(0, 0, 0);
		q->force = LR_vector<number>(0, 0, 0);
		LR_vector<number> r = q->pos.minimum_image(p->pos, *this->_config_info.box_side);

		this->_config_info.interaction->pair_interaction((*it).first, (*it).second, &r, true);

		// loop over all the 9 elements of the stress tensor
		for(int si = 0; si < 3; si++) {
			for(int sj = 0; sj < 3; sj++) {
				st[si][sj] += r[si] * q->force[sj];
			}
		}
	}

	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		for(int si = 0; si < 3; si++) {
			for(int sj = 0; sj < 3; sj++) {
				st[si][sj] += p->vel[si] * p->vel[sj];
			}
		}
	}

	stringstream res;
	if(_print_averaged_off_diagonal) {
		number avg = 0;
		for (int si = 0; si < 3; si++) {
			for (int sj = 0; sj < 3; sj++) {
				if(si != sj) avg += st[si][sj];
			}
		}
		avg /= 6;
		res << avg;
	}
	else {
		for(int si = 0; si < 3; si++) {
			for(int sj = 0; sj < 3; sj++) {
				if(si != 0 || sj != 0) res << " ";
				res << st[si][sj];
			}
		}
	}

	return res.str();
}

template class StressTensor<float>;
template class StressTensor<double>;

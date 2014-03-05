/*
 * BinaryConfiguration.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "BinaryConfiguration.h"

template<typename number>
BinaryConfiguration<number>::BinaryConfiguration() : Configuration<number>() {

}

template<typename number>
BinaryConfiguration<number>::~BinaryConfiguration() {

}

template<typename number>
std::string BinaryConfiguration<number>::_headers(llint step) {
	std::stringstream headers;
	headers.precision(15);

	OX_DEBUG("saving step: %llu", step);
	headers.write((char * )(&step), sizeof (llint));

	unsigned short rndseed[3];
	Utils::get_seed(rndseed);

	OX_DEBUG("saving seed: %hu %hu %hu", rndseed[0], rndseed[1], rndseed[2]);

	// print out the number
	headers.write((char *)rndseed, 3 * sizeof(unsigned short));

	double mybox = *this->_config_info.box_side;
	headers.write ((char * )(&mybox), sizeof (double));
	headers.write ((char * )(&mybox), sizeof (double));
	headers.write ((char * )(&mybox), sizeof (double));

	double U = (double) this->_tot_energy.get_U(step);
	double K = (double) this->_tot_energy.get_K(step);
	double tot = U+K;
	headers.write ((char * )(&tot), sizeof (double));
	headers.write ((char * )(&U), sizeof (double));
	headers.write ((char * )(&K), sizeof (double));

	return headers.str();
}

template<typename number>
std::string BinaryConfiguration<number>::_configuration(llint step) {
	std::stringstream conf;
	conf.precision(15);

	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<double> mypos(p->get_abs_pos(*this->_config_info.box_side).x, p->get_abs_pos(*this->_config_info.box_side).y, p->get_abs_pos(*this->_config_info.box_side).z);
		LR_matrix<number> oT = p->orientation.get_transpose();

		double tmpf = mypos.x;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = mypos.y;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = mypos.z;
		conf.write ((char * )(&tmpf), sizeof (double));

		tmpf = (double) oT.v1.x;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) oT.v1.y;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) oT.v1.z;
		conf.write ((char * )(&tmpf), sizeof (double));

		tmpf = (double) oT.v3.x;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) oT.v3.y;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) oT.v3.z;
		conf.write ((char * )(&tmpf), sizeof (double));

		tmpf = (double) p->vel.x;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) p->vel.y;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) p->vel.z;
		conf.write ((char * )(&tmpf), sizeof (double));

		tmpf = (double) p->L.x;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) p->L.y;
		conf.write ((char * )(&tmpf), sizeof (double));
		tmpf = (double) p->L.z;
		conf.write ((char * )(&tmpf), sizeof (double));
	}

	return conf.str();
}

template class BinaryConfiguration<float>;
template class BinaryConfiguration<double>;

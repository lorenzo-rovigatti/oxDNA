/*
 * BinaryConfiguration.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "BinaryConfiguration.h"

BinaryConfiguration::BinaryConfiguration() :
				Configuration() {

}

BinaryConfiguration::~BinaryConfiguration() {

}

std::string BinaryConfiguration::_headers(llint step) {
	std::stringstream headers;

	headers.write((char *) (&step), sizeof(llint));

	unsigned short rndseed[3];
	Utils::get_seed(rndseed);

	OX_DEBUG("Saving conf. at step %llu, rng status: %hu %hu %hu", step, rndseed[0], rndseed[1], rndseed[2]);

	// print out the number
	headers.write((char *) rndseed, 3 * sizeof(unsigned short));

	LR_vector my_box_sides(_config_info->box->box_sides().x, _config_info->box->box_sides().y, _config_info->box->box_sides().z);
	headers.write((char *) (&my_box_sides.x), sizeof(double));
	headers.write((char *) (&my_box_sides.y), sizeof(double));
	headers.write((char *) (&my_box_sides.z), sizeof(double));

	double U = (double) this->_tot_energy.get_U(step);
	double K = (double) this->_tot_energy.get_K(step);
	double tot = U + K;
	headers.write((char *) (&tot), sizeof(double));
	headers.write((char *) (&U), sizeof(double));
	headers.write((char *) (&K), sizeof(double));

	return headers.str();
}

std::string BinaryConfiguration::_configuration(llint step) {
	std::stringstream conf;

	for(int i = 0; i < _config_info->N(); i++) {
		BaseParticle *p = _config_info->particles()[i];

		LR_vector mypos(p->pos.x, p->pos.y, p->pos.z);

		LR_matrix oT = p->orientation.get_transpose();

		double tmpf = mypos.x;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = mypos.y;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = mypos.z;
		conf.write((char *) (&tmpf), sizeof(double));

		int tmpi[3];
		p->get_pos_shift((int *) tmpi);
		conf.write((char *) (&(tmpi[0])), sizeof(int));
		conf.write((char *) (&(tmpi[1])), sizeof(int));
		conf.write((char *) (&(tmpi[2])), sizeof(int));

		tmpf = (double) oT.v1.x;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) oT.v1.y;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) oT.v1.z;
		conf.write((char *) (&tmpf), sizeof(double));

		tmpf = (double) oT.v2.x;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) oT.v2.y;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) oT.v2.z;
		conf.write((char *) (&tmpf), sizeof(double));

		tmpf = (double) oT.v3.x;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) oT.v3.y;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) oT.v3.z;
		conf.write((char *) (&tmpf), sizeof(double));

		tmpf = (double) p->vel.x;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) p->vel.y;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) p->vel.z;
		conf.write((char *) (&tmpf), sizeof(double));

		tmpf = (double) p->L.x;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) p->L.y;
		conf.write((char *) (&tmpf), sizeof(double));
		tmpf = (double) p->L.z;
		conf.write((char *) (&tmpf), sizeof(double));
	}

	// renormalization part
	/* 
	 for(int i = 0; i < *_config_info->N; i++) {
	 BaseParticle *p = _config_info->particles()[i];

	 LR_matrix oT = p->orientation.get_transpose();

	 oT.v1.normalize();
	 oT.v3.normalize();

	 oT.v1 -= oT.v3 * (oT.v1 * oT.v3);
	 oT.v1.normalize();
	 oT.v2 = oT.v3.cross(oT.v1);
	 //oT.v2.normalize();

	 oT.transpone();

	 p->orientation = oT;
	 p->set_positions();
	 p->orientationT = p->orientation.get_transpose();
	 }
	 */

	return conf.str();
}

/*
 * NathanNeighs.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "NathanNeighs.h"
#include "../Utilities/Utils.h"

template<typename number>
NathanNeighs<number>::NathanNeighs() : Configuration<number>() {
	_total_bonds = 0;
	_n_crystalline = 0;
	_mode = BONDS;
}

template<typename number>
NathanNeighs<number>::~NathanNeighs() {

}

template<typename number>
void NathanNeighs<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
	getInputNumber(&my_inp, "threshold", &_threshold, 0);
	getInputNumber(&sim_inp, "NATHAN_alpha", &_patch_length, 1);
	_patch_length += 0.5;

	string my_mode;
	getInputString(&my_inp, "mode", my_mode, 1);
	if(my_mode == "mgl") _mode = MGL;
	else if(my_mode == "tot_bonds") _mode = TOT_BONDS;
	else if(my_mode == "min_theta") _mode = MIN_THETA;
	else if(my_mode == "bonds") _mode = BONDS;
	else throw oxDNAException("NathanNeighs: the only acceptable modes are mgl, tot_bonds, min_theta and bonds");
}

template<typename number>
void NathanNeighs<number>::init(ConfigInfo<number> &config_info) {
   Configuration<number>::init(config_info);
}

template<typename number>
std::string NathanNeighs<number>::_headers(llint step) {
	std::stringstream headers;

	if(_mode == MGL) {
		number mybox = *this->_config_info.box_side;
		headers << ".Box:" << mybox << "," << mybox << "," << mybox << endl;
	}

	return headers.str();
}

template<typename number>
std::string NathanNeighs<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	if(p->type != 0) return res.str();

	int n1 = 0;
	int n2 = 0;
	LR_vector<number> p_axis = p->orientationT.v3;
	// this is the number of particles which are no more than 1 + alpha far apart from p
	int n_within = 0;
	vector<BaseParticle<number> *> particles = this->_config_info.interaction->get_neighbours(p, this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);
	LR_vector<number> bonded_p1[3];
	LR_vector<number> bonded_p2[3];
	for(typename std::vector<BaseParticle<number> *>::iterator it = particles.begin(); it != particles.end(); it++) {
		BaseParticle<number> *q = *it;
		if(q->type == 0) {
			LR_vector<number> r = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
			number r_mod = r.module();
			if(r_mod < (_patch_length + 0.5)) n_within++;
			if(this->_config_info.interaction->pair_interaction(p, q) < _threshold) {
				r /= r_mod;
				if(p_axis*r > 0) {
					if(n1 < 3) {
						bonded_p1[n1] = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
						bonded_p1[n1] = bonded_p1[n1] - p_axis*(bonded_p1[n1]*p_axis);
						bonded_p1[n1].normalize();
					}
					n1++;
				}
				else {
					if(n2 < 3) {
						bonded_p2[n2] = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
						bonded_p2[n2] = bonded_p2[n2] - p_axis*(bonded_p2[n2]*p_axis);
						bonded_p2[n2].normalize();
					}
					n2++;
				}
			}
		}
	}
	if(_mode == BONDS) res << n1 << " " << n2;

	if(n1 == 3 && n2 == 3 && n_within == 6) {
		_n_crystalline++;

		number avg_min_theta = 0.;
		for(int i = 0; i < 3; i++) {
			LR_vector<number> b_pos_1 = bonded_p1[i];
			number min_theta = 1000.;
			for(int j = 0; j < 3; j++) {
				LR_vector<number> b_pos_2 = bonded_p2[j];
				number scalar_prod = b_pos_1*b_pos_2;
				if(scalar_prod < -1.) scalar_prod = -0.999999999999;
				else if(scalar_prod > 1.) scalar_prod = 0.9999999999;
				number theta = acos(scalar_prod);
				if(theta < min_theta) min_theta = theta;
			}
			avg_min_theta += min_theta;
		}
		avg_min_theta /= 3.;
		_avg_min_theta += avg_min_theta;

		if(_mode == MGL) {
			LR_vector<number> p1 = p_axis*_patch_length;
			LR_vector<number> p2 = -p_axis*_patch_length;
			string str = Utils::sformat("%lf %lf %lf @ 0.5 C[red] M %lf %lf %lf %lf C[blue] %lf %lf %lf %lf C[blue]", p->pos.x, p->pos.y,p->pos.z, p1.x, p1.y, p1.z, 0.7, p2.x, p2.y, p2.z, 0.7);
			res << str;
		}
		if(_mode == MIN_THETA) res << avg_min_theta;
	}
	_total_bonds += n1 + n2;

	return res.str();
}

template<typename number>
std::string NathanNeighs<number>::_configuration(llint step) {
	_total_bonds = 0;
	_n_crystalline = 0;
	_avg_min_theta = 0.;
	std::string tot_conf = Configuration<number>::_configuration(step);

	if(_mode == TOT_BONDS) {
		double avg_theta = (_n_crystalline > 0) ? _avg_min_theta / (number)_n_crystalline: 0.;
		return Utils::sformat("%d %lf %d", _total_bonds, avg_theta, _n_crystalline);
	}
	else return tot_conf;
}

template class NathanNeighs<float>;
template class NathanNeighs<double>;

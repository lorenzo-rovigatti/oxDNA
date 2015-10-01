/*
 * NathanNeighs.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include "NathanNeighs.h"
#include "../Utilities/Utils.h"

#include <sstream>
#include <cmath>

#include <gsl/gsl_sf.h>

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
	else if(my_mode == "qs") _mode = QS;
	else if(my_mode == "qs_avg") _mode = QS_AVG;
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
void NathanNeighs<number>::_compute_qn(NateParticle<number> &np, int n) {
	vector<complex<number> > qm(2*n + 1);
	number norm = 0.;

	for(int m = -n; m <= n; m++) {
		int mm = m + n;
		int am = abs(m);
		qm[mm] = 0.;

		typename vector<LR_vector<number> >::iterator it;
		for(it = np.vs.begin(); it != np.vs.end(); it++) {
			LR_vector<number> v = *it;
			number Plm = gsl_sf_legendre_sphPlm(n, am, v[2]);
			number phi = atan2(v[1], v[0]);
			if(phi < 0.) phi += 2*M_PI;

			complex<number> qtemp(Plm*cos(am*phi), Plm*sin(am*phi));

			if(m < 0) {
				qtemp = conj(qtemp);
				qtemp *= pow(-1, am);
			}

			qm[mm] += qtemp;
		}

		qm[mm] /= (number) np.vs.size();
		norm += std::norm(qm[mm]);
	}

	// normalise the vector
	norm = sqrt(norm);
	for(int m = -n; m <= n; m++) {
		int mm = n + m;
		qm[mm] /= norm;
	}

	if(n == 4) np.q4 = qm;
	else np.q6 = qm;
}

template<typename number>
std::string NathanNeighs<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	if(p->type != 0) return res.str();
	
	NateParticle<number> &np = _nate_particles[p->index];
	np.is_crystalline = false;
	np.p = p;

	int n1 = 0;
	int n2 = 0;
	LR_vector<number> p_axis = p->orientationT.v3;
	// this is the number of particles which are no more than 1 + alpha far apart from p
	int n_within = 0;
	vector<BaseParticle<number> *> particles = this->_config_info.lists->get_all_neighbours(p);
	for(typename std::vector<BaseParticle<number> *>::iterator it = particles.begin(); it != particles.end(); it++) {
		BaseParticle<number> *q = *it;
		if(q->type == 0) {
			LR_vector<number> r = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
			number r_mod = r.module();
			if(r_mod < (_patch_length + 0.5)) n_within++;
			if(this->_config_info.interaction->pair_interaction(p, q) < _threshold) {
				np.neighs.push_back(&_nate_particles[q->index]);
				r /= r_mod;
				if(p_axis*r > 0) {
					if(n1 < 3) {
						LR_vector<number> v = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
						v.normalize();
						np.vs.push_back(v);
						v = v - p_axis*(v*p_axis);
						v.normalize();
						np.v1.push_back(v);
					}
					n1++;
				}
				else {
					if(n2 < 3) {
						LR_vector<number> v = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
						v.normalize();
						np.vs.push_back(v);
						v = v - p_axis*(v*p_axis);
						v.normalize();
						np.v2.push_back(v);
					}
					n2++;
				}
			}
		}
	}
	if(_mode == BONDS) res << n1 << " " << n2;

	_compute_qn(np, 4);
	_compute_qn(np, 6);

	if(n1 == 3 && n2 == 3 && n_within == 6) {
		np.is_crystalline = true;
		_n_crystalline++;

		number avg_min_theta = 0.;
		for(int i = 0; i < 3; i++) {
			LR_vector<number> b_pos_1 = np.v1[i];
			number min_theta = 1000.;
			for(int j = 0; j < 3; j++) {
				LR_vector<number> b_pos_2 = np.v2[j];
				number scalar_prod = b_pos_1*b_pos_2;
				if(scalar_prod < -1.) scalar_prod = -0.999999999999;
				else if(scalar_prod > 1.) scalar_prod = 0.99999999999;
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
	_nate_particles.clear();
	_nate_particles.resize(*this->_config_info.N);

	_total_bonds = 0;
	_n_crystalline = 0;
	_avg_min_theta = 0.;
	std::string tot_conf = Configuration<number>::_configuration(step);

	if(_mode == TOT_BONDS) {
		double avg_theta = (_n_crystalline > 0) ? _avg_min_theta / (number)_n_crystalline: 0.;
		return Utils::sformat("%d %lf %d", _total_bonds, avg_theta, _n_crystalline);
	}
	else if(_mode == QS) {
		stringstream ss;
		bool first = true;
		for(int i = 0; i < *this->_config_info.N; i++) {
			NateParticle<number> *np = &_nate_particles[i];
			if(np->is_crystalline) {
				typename vector<NateParticle<number> *>::iterator it;
				for(it = np->neighs.begin(); it != np->neighs.end(); it++) {
					NateParticle<number> *nq = *it;
					// q4
					complex<number> qiqj;
					for(uint di = 0; di < np->q4.size(); di++) qiqj += np->q4[di]*conj(nq->q4[di]);
					if(!first) ss << endl;
					else first = false;
					ss << qiqj.real();
					// q6
					qiqj = complex<number>();
					for(uint di = 0; di < np->q6.size(); di++) qiqj += np->q6[di]*conj(nq->q6[di]);
					ss << " " << qiqj.real();
				}
			}
		}

		return ss.str();
	}
	else if(_mode == QS_AVG) {
		stringstream ss;
		ss << "# " << _n_crystalline;
		for(int i = 0; i < *this->_config_info.N; i++) {
			NateParticle<number> *np = &_nate_particles[i];
			if(np->is_crystalline) {
				vector<complex<number> > q4(np->q4);
				vector<complex<number> > q6(np->q6);
				typename vector<NateParticle<number> *>::iterator it;
				for(it = np->neighs.begin(); it != np->neighs.end(); it++) {
					for(int m = -4; m <= 4; m++) {
						int mm = 4 + m;
						q4[mm] += (*it)->q4[mm];
					}
					for(int m = -6; m <= 6; m++) {
						int mm = 6 + m;
						q6[mm] += (*it)->q6[mm];
					}
				}

				number avg_q4 = 0.;
				for(int m = -4; m <= 4; m++) {
					int mm = 4 + m;
					avg_q4 += std::norm(q4[mm]);
				}
				avg_q4 *= 4*M_PI/(9.*SQR(np->neighs.size()));

				number avg_q6 = 0.;
				for(int m = -6; m <= 6; m++) {
					int mm = 6 + m;
					avg_q6 += std::norm(q6[mm]);
				}
				avg_q6 *= 4*M_PI/(13.*SQR(np->neighs.size()));

				ss << endl << np->p->index << " " << avg_q4 << " " << avg_q6;
			}
		}

		return ss.str();
	}
	else return tot_conf;
}

template class NathanNeighs<float>;
template class NathanNeighs<double>;

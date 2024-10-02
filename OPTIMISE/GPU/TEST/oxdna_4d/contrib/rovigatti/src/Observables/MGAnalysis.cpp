/*
 * ElasticConstantTensor.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include "MGAnalysis.h"

#include "Utilities/Utils.h"

#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"
#include "diagonalise_3x3.h"

MGAnalysis::MGAnalysis() :
				_exponent(6),
				_alpha(0.),
				_volume_only(true),
				_rg_only(false) {
	_two_microgels = false;
}

MGAnalysis::~MGAnalysis() {

}

void MGAnalysis::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	std::string inter;
	getInputString(&sim_inp, "interaction_type", inter, 1);
	if(inter != "MGInteraction") throw oxDNAException("ElasticConstantTensor is not compatible with the interaction '%s'", inter.c_str());

	getInputInt(&sim_inp, "MG_n", &_exponent, 1);

	number rcut;
	getInputNumber(&sim_inp, "MG_rcut", &rcut, 1);
	_sqr_rcut = SQR(rcut);

	getInputNumber(&sim_inp, "MG_alpha", &_alpha, 0);

	number rfene;
	getInputNumber(&sim_inp, "MG_rfene", &rfene, 1);
	_sqr_rfene = SQR(rfene);

	getInputBool(&my_inp, "volume_only", &_volume_only, 0);
	getInputBool(&my_inp, "rg_only", &_rg_only, 0);
	getInputBool(&my_inp, "two_microgels", &_two_microgels, 0);

	if(_volume_only) {
		_volume_threshold = 0.;
		getInputNumber(&my_inp, "volume_threshold", &_volume_threshold, 0);
		_volume_threshold_sqr = SQR(_volume_threshold);
	}

	if(_volume_only && _rg_only) {
		throw oxDNAException("MGAnalysis: volume_only and rg_only are incompatible");
	}

	if(_two_microgels) {
		OX_LOG(Logger::LOG_INFO, "MGAnalysis: Assuming configurations with two identical microgels");
	}
}

void MGAnalysis::init() {
	number rep_rcut = pow(2., 1. / _exponent);
	_sqr_rep_rcut = SQR(rep_rcut);

	_gamma = M_PI / (2.25 - pow(2., 1. / 3.));
	_beta = 2 * M_PI - 2.25 * _gamma;
}

std::pair<number, number> MGAnalysis::_lame_coefficients() {
	std::vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();

	number lambda = 0.;
	number mu = 0.;
	for(auto &pair: pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;
		LR_vector r = _config_info->box->min_image(p->pos, q->pos);
		number r_sqr = r.norm();
		number r_mod = sqrt(r_sqr);

		number first_der = 0.;
		number second_der = 0.;

		if(p->is_bonded(q)) {
			first_der += 30 * _sqr_rfene * r_mod / (_sqr_rfene - r_sqr);
			second_der += 30 * _sqr_rfene / (_sqr_rfene - r_sqr) * (1. + 2 * r_sqr / (_sqr_rfene - r_sqr));
		}

		if(r_sqr < _sqr_rcut) {
			if(r_sqr < _sqr_rep_rcut) {
				first_der += 4. * (6. * pow(r_mod, -7.) - 12. * pow(r_mod, -13.));
				second_der += 4. * (12. * 13. * pow(r_mod, -14.) - 6. * 7. * pow(r_mod, -8.));
			}
			else {
				number sin_part = sin(_gamma * r_sqr + _beta);
				number cos_part = cos(_gamma * r_sqr + _beta);
				first_der += -_alpha * _gamma * r_mod * sin_part;
				second_der += -(_alpha * _gamma * sin_part + 2 * _alpha * SQR(_gamma) * r_sqr * cos_part);
			}
		}

		number factor = second_der - first_der / r_mod;
		for(int d1 = 0; d1 < 3; d1++) {
			for(int d2 = 0; d2 < 3; d2++) {
				if(d1 != d2) {
					number contrib = factor * SQR(r[d1]) * SQR(r[d2]) / r_sqr;
					lambda += contrib;
					mu += contrib;
				}
			}
		}
	}

	lambda /= 6.;
	mu /= 6.;
	mu += 2 * _config_info->temperature() * _config_info->N();

	return std::pair<number, number>(lambda, mu);
}

number MGAnalysis::_volume() {
	int N = _config_info->N();
	LR_vector com = _com(0, N);

	std::vector<qh_vertex_t> vertices(N);
	int curr_idx = 0;
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		LR_vector p_pos = _config_info->box->get_abs_pos(p) - com;

		if(_volume_threshold == 0. || p_pos.norm() < _volume_threshold_sqr) {
			vertices[curr_idx].x = p_pos.x;
			vertices[curr_idx].y = p_pos.y;
			vertices[curr_idx].z = p_pos.z;
			curr_idx++;
		}
	}

	qh_mesh_t mesh = qh_quickhull3d(vertices.data(), curr_idx);

	number volume = 0.;
	for(int i = 0, j = 0; i < (int) mesh.nindices; i += 3, j++) {
		LR_vector p1(mesh.vertices[mesh.indices[i + 0]].x, mesh.vertices[mesh.indices[i + 0]].y, mesh.vertices[mesh.indices[i + 0]].z);
		LR_vector p2(mesh.vertices[mesh.indices[i + 1]].x, mesh.vertices[mesh.indices[i + 1]].y, mesh.vertices[mesh.indices[i + 1]].z);
		LR_vector p3(mesh.vertices[mesh.indices[i + 2]].x, mesh.vertices[mesh.indices[i + 2]].y, mesh.vertices[mesh.indices[i + 2]].z);

		volume += (p1 * (p2.cross(p3))) / 6.;
	}

	qh_free_mesh(mesh);

	return volume;
}

LR_vector MGAnalysis::_com(int from_idx, int to_idx) {
	LR_vector com;
	int N = to_idx - from_idx;
	for(int i = from_idx; i < to_idx; i++) {
		BaseParticle *p = _config_info->particles()[i];
		com += _config_info->box->get_abs_pos(p);
	}
	return com / N;
}

std::vector<number> MGAnalysis::_rg_eigenvalues(int from_idx, int to_idx) {
	std::vector<number> res;
	LR_vector com = _com(from_idx, to_idx);
	int N = to_idx - from_idx;

	double IM[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	for(int i = from_idx; i < to_idx; i++) {
		BaseParticle *p = _config_info->particles()[i];
		LR_vector i_pos = _config_info->box->get_abs_pos(p) - com;

		IM[0][0] += SQR(i_pos[1]) + SQR(i_pos[2]);
		IM[0][1] += -i_pos[0] * i_pos[1];
		IM[0][2] += -i_pos[0] * i_pos[2];

		IM[1][1] += SQR(i_pos[0]) + SQR(i_pos[2]);
		IM[1][2] += -i_pos[1] * i_pos[2];

		IM[2][2] += SQR(i_pos[0]) + SQR(i_pos[1]);
	}
	IM[1][0] = IM[0][1];
	IM[2][0] = IM[0][2];
	IM[2][1] = IM[1][2];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			IM[i][j] /= N;

	double EV[3][3];
	double val[3];
	eigen_decomposition(IM, EV, val);

	res.push_back(sqrt(val[0]));
	res.push_back(sqrt(val[1]));
	res.push_back(sqrt(val[2]));

	LR_vector EVs[3] = { LR_vector(EV[0][0], EV[0][1], EV[0][2]), LR_vector(EV[1][0], EV[1][1], EV[1][2]), LR_vector(EV[2][0], EV[2][1], EV[2][2]) };
	EVs[0].normalize();
	EVs[1].normalize();
	EVs[2].normalize();

	LR_vector max_along_EVs(-1.e6, -1.e6, -1.e6);
	LR_vector min_along_EVs(1.e6, 1.e6, 1.e6);

	for(int i = from_idx; i < to_idx; i++) {
		BaseParticle *p = _config_info->particles()[i];
		LR_vector p_pos = _config_info->box->get_abs_pos(p) - com;

		for(int d = 0; d < 3; d++) {
			number abs = p_pos * EVs[d];
			if(abs > max_along_EVs[d]) max_along_EVs[d] = abs;
			else if(abs < min_along_EVs[d]) min_along_EVs[d] = abs;
		}
	}

	res.push_back(max_along_EVs[0] - min_along_EVs[0]);
	res.push_back(max_along_EVs[1] - min_along_EVs[1]);
	res.push_back(max_along_EVs[2] - min_along_EVs[2]);

	return res;
}

std::string MGAnalysis::get_output_string(llint curr_step) {
	std::string to_ret;

	if(_volume_only) {
		double volume = _volume();
		to_ret = Utils::sformat("%lf", volume);
	}
	else if(_rg_only) {
		int N = _config_info->N();
		if(!_two_microgels) {
			auto eigenvalues = _rg_eigenvalues(0, N);
			to_ret = Utils::sformat("%lf %lf %lf %lf %lf %lf", eigenvalues[0], eigenvalues[1], eigenvalues[2], eigenvalues[3], eigenvalues[4], eigenvalues[5]);
		}
		else {
			int N_half = N / 2;
			auto eigenvalues = _rg_eigenvalues(0, N_half);
			to_ret = Utils::sformat("%lf %lf %lf", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
			eigenvalues = _rg_eigenvalues(N_half, N);
			to_ret += Utils::sformat(" %lf %lf %lf", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
		}
	}
	else {
		double volume = _volume();
		std::pair<number, number> lame = _lame_coefficients();
		number lambda = lame.first / volume;
		number mu = lame.second / volume;
		to_ret = Utils::sformat("%lf %lf", lambda, mu);
	}

	return to_ret;
}

/*
 * PatchyInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "PatchyInteraction.h"
#include "../Particles/PatchyParticle.h"
#include "../Utilities/Utils.h"

template <typename number>
PatchyInteraction<number>::PatchyInteraction() : BaseInteraction<number, PatchyInteraction<number> >(), _N_patches(0), _N_patches_B(-1), _N_A(0), _N_B(0), _is_binary(false) {
	this->_int_map[PATCHY] = &PatchyInteraction<number>::_patchy_interaction;

	for(int i = 0; i < 3; i++) {
		_sigma[i] = 1.;
		_sqr_sigma[i] = 1.;
		_epsilon[i] = 1.;
	}

	_patch_alpha = 0.12;
}

template <typename number>
PatchyInteraction<number>::~PatchyInteraction() {

}

template<typename number>
void PatchyInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputInt(&inp, "PATCHY_N", &_N_patches, 1);
	if(getInputInt(&inp, "PATCHY_N_B", &_N_patches_B, 0) == KEY_FOUND) _is_binary = true;

	if(_is_binary) {
		getInputNumber<number>(&inp, "PATCHY_sigma_AA", _sigma, 0);
		getInputNumber<number>(&inp, "PATCHY_sigma_BB", _sigma + 2, 0);
		if(getInputNumber<number>(&inp, "PATCHY_sigma_AB", _sigma + 1, 0) == KEY_NOT_FOUND) {
			_sigma[1] = (_sigma[0] + _sigma[2])*0.5;
		}

		getInputNumber<number>(&inp, "PATCHY_epsilon_AA", _epsilon, 0);
		getInputNumber<number>(&inp, "PATCHY_epsilon_BB", _epsilon + 2, 0);
		if(getInputNumber<number>(&inp, "PATCHY_epsilon_AB", _epsilon + 1, 0) == KEY_NOT_FOUND) {
			_epsilon[1] = sqrt(_epsilon[0]*_epsilon[2]);
		}
	}

	getInputNumber<number>(&inp, "PATCHY_alpha", &_patch_alpha, 0);
}

template<typename number>
void PatchyInteraction<number>::init() {
	number patch_rcut = _patch_alpha*1.5;
	_sqr_patch_rcut = SQR(patch_rcut);
	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);
	number r8b10 = powf(patch_rcut, (number) 8.f) / _patch_pow_alpha;
	this->_rcut = 0;
	for(int i = 0; i < 3; i++) {
		number rcut = _sigma[i]*1.05 + patch_rcut;
		if(rcut > this->_rcut) this->_rcut = rcut;
		_sqr_tot_rcut[i] = SQR(rcut);
		_sqr_sigma[i] = SQR(_sigma[i]);
		_patch_E_cut[i] = -1.001f*_epsilon[i]*expf(-(number)0.5f*r8b10*_sqr_patch_rcut);
		_E_cut[i] = powf((number) _sigma[i]/rcut, PATCHY_POWER);
	}

	this->_sqr_rcut = SQR(this->_rcut);

	if(_is_binary) OX_LOG(Logger::LOG_INFO, "Simulating a binary patchy system with diameters AA %lf, BB %lf, AB %lf (N patch: %d %d, rcut: %lf %lf %lf)", _sigma[0], _sigma[2], _sigma[1], _N_patches, _N_patches_B, sqrt(_sqr_tot_rcut[0]), sqrt(_sqr_tot_rcut[1]), sqrt(_sqr_tot_rcut[2]));
	else OX_LOG(Logger::LOG_INFO, "Simulating a pure patchy system (N patch: %d, rcut: %lf, patch_alpha: %lf\n", _N_patches, this->_rcut);
}

template<typename number>
void PatchyInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		if(i < _N_A) particles[i] = new PatchyParticle<number>(_N_patches, P_A, _sigma[2*P_A]);
		else particles[i] = new PatchyParticle<number>(_N_patches_B, P_B, _sigma[2*P_B]);
	}
}

template<typename number>
number PatchyInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number PatchyInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number PatchyInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _patchy_interaction(p, q, r, update_forces);
}

template<typename number>
void PatchyInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_B);
	if(_N_B > 0) if(_N_patches_B == -1) throw oxDNAException("Number of patches of species B not specified");
	_N_A = N - _N_B;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = (i < _N_A) ? P_A : P_B;
	   particles[i]->btype = (i < _N_A) ? P_A : P_B;
	   particles[i]->strand_id = i;
	}
}

template<typename number>
void PatchyInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class PatchyInteraction<float>;
template class PatchyInteraction<double>;

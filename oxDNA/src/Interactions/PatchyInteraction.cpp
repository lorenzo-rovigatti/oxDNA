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
PatchyInteraction<number>::PatchyInteraction() : BaseInteraction<number, PatchyInteraction<number> >(), _N_patches_B(-1), _N_B(0), _is_binary(false) {
	this->_int_map[PATCHY] = &PatchyInteraction<number>::_patchy_interaction;
}

template <typename number>
PatchyInteraction<number>::~PatchyInteraction() {

}

template<typename number>
void PatchyInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputInt(&inp, "PATCHY_N", &_N_patches, 1);
	if(getInputInt(&inp, "PATCHY_N_B", &_N_patches_B, 0) == KEY_FOUND) _is_binary = true;

	float tmp = 1.2;
	getInputFloat(&inp, "PATCHY_rcut", &tmp, 0);
	this->_rcut = (number) tmp;

	tmp = 0.12;
	getInputFloat(&inp, "PATCHY_alpha", &tmp, 0);
	_patch_alpha = (number) tmp;
}

template<typename number>
void PatchyInteraction<number>::init() {
	_E_cut = powf((number) this->_rcut, -PATCHY_POWER);
	this->_sqr_rcut = SQR(this->_rcut);

	_patch_pow_alpha = powf(_patch_alpha, (number) 10.f);

	number r8b10 = powf(PATCHY_CUTOFF, (number) 8.f) / _patch_pow_alpha;
	_patch_E_cut = -1.001f * expf(-(number)0.5f * r8b10 * SQR(PATCHY_CUTOFF));
}

template<typename number>
void PatchyInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		int i_patches = (i < _N_A) ? _N_patches : _N_patches_B;
		particles[i] = new PatchyParticle<number>(i_patches);
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

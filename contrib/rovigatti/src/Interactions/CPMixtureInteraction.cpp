/*
 * CPMixtureInteraction.cpp
 *
 *  Created on: 05/may/2016
 *      Author: lorenzo
 */

#include "CPMixtureInteraction.h"

template<typename number>
CPMixtureInteraction<number>::CPMixtureInteraction() : BaseInteraction<number, CPMixtureInteraction<number> >() {
	this->_int_map[CP] = &CPMixtureInteraction<number>::pair_interaction_nonbonded;
	_sigma[0] = _sigma[1] = _sigma[2] = 1.;
	_epsilon[0] = _epsilon[1] = _epsilon[2] = 1.;
	_n[0] = _n[1] = _n[2] = 6;
	_N_A = _N_B = 0;
}

template<typename number>
CPMixtureInteraction<number>::~CPMixtureInteraction() {

}

template<typename number>
void CPMixtureInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	 for(int i = 0; i < 3; i++) {
			char key[256];
			sprintf(key, "CP_sigma[%d]", i);
			getInputNumber(&inp, key, &_sigma[i], 1);
			sprintf(key, "CP_rcut[%d]", i);
			getInputNumber(&inp, key, &_CP_rcut[i], 1);
			sprintf(key, "CP_epsilon[%d]", i);
			getInputNumber(&inp, key, &_epsilon[i], 1);
			sprintf(key, "CP_type[%d]", i);
			string int_type;
			getInputString(&inp, key, int_type, 1);

			if(int_type == "power_law") {
				_CP_int_type[i] = POWER_LAW;
				sprintf(key, "CP_n[%d]", i);
				getInputInt(&inp, key, _n + i, 1);
				// no odd exponents allowed
				if((_n[0] % 2) == 1) throw oxDNAException("%s should be an even integer", key);
			}
			else if(int_type == "gaussian") _CP_int_type[i] = GAUSSIAN;
			else if(int_type == "louis") _CP_int_type[i] = LOUIS;
			else throw oxDNAException("CP interaction type '%s' not supported", int_type.c_str());
	}
}

template<typename number>
void CPMixtureInteraction<number>::init() {
	this->_rcut = 0;
	for(int i = 0; i < 3; i++) {
		_CP_rcut[i] *= _sigma[i];
		if(_CP_rcut[i] > this->_rcut) this->_rcut = _CP_rcut[i];
		_sqr_CP_rcut[i] = SQR(_CP_rcut[i]);
		_sqr_sigma[i] = SQR(_sigma[i]);

		OX_LOG(Logger::LOG_INFO, "CPMixture, interaction n. %d, type = %d, r_cut = %lf", i, _CP_int_type[i], _CP_rcut[i]);

		switch(_CP_int_type[i]) {
		case POWER_LAW:
			_E_cut[i] = _epsilon[i]*pow(_sigma[i]/_CP_rcut[i], _n[i]);
			break;
		case GAUSSIAN:
			 _E_cut[i] = _epsilon[i]*exp(-4*_sqr_CP_rcut[i]/_sqr_sigma[i]);
			break;
		case LOUIS: {
			_A[i] = 401.;
			_C[i] = -159.;
			number exp_factor = 2./_sigma[i]/0.8;
			_B[i] = -3.686*exp_factor;
			_D[i] = -3.28*pow(exp_factor, 1.5);

			_E_cut[i] = _epsilon[i]*(_A[i]*exp(_B[i]*_CP_rcut[i]) + _C[i]*exp(_D[i]*pow(_CP_rcut[i], 1.5)));
			break;
		}
		default:
			throw oxDNAException("I shouldn't be here...");
			break;
		}
	}

	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void CPMixtureInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
void CPMixtureInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_B);
	_N_A = N - _N_B;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = particles[i]->strand_id = i;
	   particles[i]->type = particles[i]->btype = (i < _N_A) ? P_A : P_B;
	}
}

template<typename number>
number CPMixtureInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number CPMixtureInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number CPMixtureInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	number sqr_r = r->norm();
	number energy = 0.;

	int type = p->type + q->type;
	if(sqr_r < _sqr_CP_rcut[type]) {
		energy = _E_cut[type];
		number force_mod = 0.;

		switch(_CP_int_type[type]) {
		case POWER_LAW:
			energy += _epsilon[type]*pow(_sqr_sigma[type]/sqr_r, _n[type]/2);
			force_mod = _n[type]*energy/sqr_r;
			break;
		case GAUSSIAN:
			energy  += _epsilon[type]*exp(-4*sqr_r/_sqr_sigma[type]);
			force_mod = 8*energy/_sqr_sigma[type];
			break;
		case LOUIS: {
			number mod_r = sqrt(sqr_r);
			number sqrt_r = sqrt(mod_r);
			number exp1 = _epsilon[type]*_A[type]*exp(_B[type]*mod_r);
			number exp2 = _epsilon[type]*_C[type]*exp(_D[type]*mod_r*sqrt_r);
			energy += exp1 + exp2;
			force_mod = -(_B[type]*exp1/mod_r + 1.5*_D[type]*exp2/sqrt_r);
			break;
		}
		default:
			throw oxDNAException("I shouldn't be here...");
			break;
		}

		if(update_forces) {
			LR_vector<number> tot_force = computed_r*force_mod;
			p->force -= tot_force;
			q->force += tot_force;
		}
	}

	return energy;
}

template<typename number>
void CPMixtureInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

extern "C" CPMixtureInteraction<float> *make_CPMixtureInteraction_float() {
	return new CPMixtureInteraction<float>();
}

extern "C" CPMixtureInteraction<double> *make_CPMixtureInteraction_double() {
	return new CPMixtureInteraction<double>();
}

template class CPMixtureInteraction<float>;
template class CPMixtureInteraction<double>;

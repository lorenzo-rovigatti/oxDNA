/*
 * JqordanInteraction.cpp
 *
 *  Created on: 19/feb/2016
 *      Author: flavio
 */

#include "JordanInteraction.h"
#include "../Particles/JordanParticle.h"
#include "../Utilities/Utils.h"

template <typename number>
JordanInteraction<number>::JordanInteraction() : BaseInteraction<number, JordanInteraction<number> >() {
	this->_int_map[JORDAN] = &JordanInteraction<number>::_jordan_interaction;

	_s = (number) 0.3; // patch width
	_m = 6;            // LJ exponent
    _phi = M_PI / 6.;  // angle below the equator
	_int_k = 0.f;      // stiffness of the internal spring
	_SB_s = 1.;       
	_n_patches = -1;
	
	_my_N = -1;
	_my_N3 = -1;
}

template <typename number>
JordanInteraction<number>::~JordanInteraction() {

}

template<typename number>
void JordanInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	float tmp = 2.5;
	getInputFloat(&inp, "JORDAN_rcut", &tmp, 0);
	getInputInt(&inp, "JORDAN_N_patches", &_n_patches, 1);
	this->_rcut = (number) tmp;
	this->_sqr_rcut = this->_rcut * this->_rcut;

	getInputNumber(&inp, "JORDAN_phi", &_phi, 0);
	getInputNumber(&inp, "JORDAN_s", &_s, 0);
	getInputNumber(&inp, "JORDAN_int_k", &_int_k, 0);
	getInputNumber(&inp, "JORDAN_SB_s", &_SB_s, 0);
	getInputInt(&inp, "JORDAN_m", &_m, 0);
	//getInputString(&inp, "JORDAN_load_patches_from", patches_file_in, 0);
	//getInputString(&inp, "JORDAN_save_patches_to", patches_file_out, 1);
}

template<typename number>
void JordanInteraction<number>::init() {
	OX_LOG(Logger::LOG_INFO,"(JordanInteraction.cpp) Running Jordan interaction with s=%g, rcut=%g, m=%d, phi=%g, int_k=%g", _s, this->_rcut, _m, _phi, _int_k, _my_N3);
}

template<typename number>
void JordanInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	OX_LOG(Logger::LOG_INFO,"(JordanInteraction.cpp) Allocating %d particles with 3 patches and %d with 4 patches", _my_N3, N - _my_N3);
	for(int i = 0; i < N; i++) {
		if (i > _my_N3 && _my_N3 >= 0) particles[i] = new JordanParticle<number>(4, _phi, _int_k);
		else particles[i] = new JordanParticle<number>(3, _phi, _int_k);
	}
}

template<typename number>
number JordanInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number JordanInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number JordanInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _jordan_interaction(p, q, r, update_forces);
}

template <typename number>
number JordanInteraction<number>::_jordan_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (update_forces) throw oxDNAException ("JordanInteraction cannot compute forces. Use MC/MC2/VMMC to simulate.");

	number rnorm = r->norm();

	// if outside of rcut, we return here
	if(rnorm > this->_sqr_rcut) return (number) 0.f;

	number rm = sqrt(rnorm);

	// otherwise, we keep going;
	number r_to_m = pow (rm, - (number)_m); // this returns (1./|r|)^m
	
	number energy = 4.f * (r_to_m * r_to_m - r_to_m);

	if (rm <= 1.f) return energy; // we do not consider the patches if r < 1

	// here we look for the closest patch on particle p;
	int imax1 = -1, imax2 = -1;
	number max1 = -2.f;
	for (int i = 0; i < p->N_int_centers; i ++) {
		number tmp = ((p->int_centers[i] * (*r)) / (0.5 * rm));
		if (tmp > max1) {
			max1 = tmp;
			imax1 = i;
		}
	}
	number max2 = -2.f;
	for (int i = 0; i < q->N_int_centers; i ++) {
		number tmp = (-(q->int_centers[i] * (*r)) / (0.5 * rm));
		if (tmp > max2) {
			max2 = tmp;
			imax2 = i;
		}
	}

	// now max1 and max2 are ready to be acos'd
	if (fabs(max1) >= 1.) max1 = copysign(1., max1);
	if (fabs(max2) >= 1.) max2 = copysign(1., max2);

	number theta1 = acos(max1);
	number theta2 = acos(max2);

	number fact = exp ((- theta1 * theta1 - theta2 * theta2) / (2.f * _s * _s));
	
	// simmetry breaking term
	// we get a random vetor and remove its projection onto r
	LR_vector<number> pprime = p->int_centers[imax1] - (*r) * ((p->int_centers[imax1] * (*r)) / (rm * rm));
	LR_vector<number> qprime = q->int_centers[imax2] - (*r) * ((q->int_centers[imax2] * (*r)) / (rm * rm));
	
	/* // check that pprime and qprime have no residual component on r
	number chk1 = pprime * (*r);
	number chk2 = qprime * (*r);
	if (fabs(chk1) > 1.e-12) printf ("Non ci siamo... %g\n", chk1);
	if (fabs(chk2) > 1.e-12) printf ("Non ci siamo proprio... %g\n", chk2);
	*/

	number Gamma_arg = (pprime * qprime) * (1.f / (pprime.module() * qprime.module()));
	if (fabs(Gamma_arg) > 1.f) Gamma_arg = copysign(1.f, Gamma_arg);
	number Gamma = acos(Gamma_arg);

	number fact_SB = exp (- Gamma * Gamma / (2.f * _SB_s * _SB_s));  // could be optimized joining with previous exponential

	return energy * fact * fact_SB;

}

template<typename number>
void JordanInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("(JordanInteraction.cpp) Can't read topology file '%s'. Aborting", this->_topology_filename);
	std::string line;
	std::getline (topology, line);
	int check = sscanf(line.c_str(), "%d %d\n", &_my_N, &_my_N3);
	if (check != 2) throw oxDNAException("(JordanInteraction.cpp) Wrong content in topology file '%s'. Expecting <N> and <N3>", this->_topology_filename);
	topology.close();

	*N_strands = N;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = 0;
	   particles[i]->btype = 1;
	   particles[i]->strand_id = i;
	}
}

template<typename number>
void JordanInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class JordanInteraction<float>;
template class JordanInteraction<double>;


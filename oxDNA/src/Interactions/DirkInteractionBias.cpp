/*
 * DirkInteractionBias.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "DirkInteractionBias.h"

template<typename number>
DirkInteractionBias<number>::DirkInteractionBias() : BaseInteraction<number, DirkInteractionBias<number> >() {
	this->_int_map[HARD] = &DirkInteractionBias<number>::_hard;
	this->_int_map[DIPOLAR] = &DirkInteractionBias<number>::_dipolar;
	this->_int_map[CHAIN] = &DirkInteractionBias<number>::_chain; 
	this->_int_map[HARM] = &DirkInteractionBias<number>::_bias;

	_K = 0.;
	_r0 = 1.;
	_chain_stiff = -1.;
	_chain_r0 = 1.;
}

template<typename number>
DirkInteractionBias<number>::~DirkInteractionBias() {

}

template<typename number>
void DirkInteractionBias<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512) && strncmp(tmps, "MC2", 512)) throw oxDNAException ("Cannot run Dirk with MD");
	
	// length
	float tmpf;
	getInputFloat (&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if (_length < 0.) throw oxDNAException ("Cannot run Dirk interactions with negative lenght");
	
	getInputFloat (&inp, "DHS_radius", &tmpf, 1);
	_DHS_radius = (number) tmpf;
	if (_DHS_radius < 0.) throw oxDNAException ("Cannot run Dirs's interaction with a negative DHS_radius");
	
	getInputFloat (&inp, "DHS_eps", &tmpf, 1);
	_DHS_eps = (number) tmpf;
	if (_DHS_eps < 0.) throw oxDNAException ("Cannot run Dirks's interaction with a negative DHS_eps");
	
	getInputFloat (&inp, "DHS_rcut", &tmpf, 1);
	_DHS_rcut = (number) tmpf;
	if (_DHS_rcut < 0.) throw oxDNAException ("Cannot run Dirks's interaction with a negative DHS_rcut");
	_DHS_sqr_rcut = _DHS_rcut * _DHS_rcut;
	
	_DHS_rf_fact = (_DHS_eps - 1.) / (2. * _DHS_eps + 1.) / (_DHS_rcut * _DHS_rcut * _DHS_rcut);
	
	_hard_rcut = 1.001 * (2. * _DHS_radius + _length);
	_hard_sqr_rcut = _hard_rcut * _hard_rcut;
	if (_hard_rcut > _DHS_rcut) this->_rcut = _hard_rcut;
	else this->_rcut = _DHS_rcut; 
	
	if (_DHS_rcut < _length) throw oxDNAException ("can't do...\n");

	getInputFloat (&inp, "K", &tmpf, 1);
	_K = (number) tmpf;
	if (_K < 0.) throw oxDNAException ("(DirkInteractionBias.cpp) Can't run with K < 0.; Aborting");
	
	getInputFloat (&inp, "r0", &tmpf, 1);
	_r0 = (number) tmpf;
	if (_r0 < 0.) throw oxDNAException ("(DirkInteractionBias.cpp) Can't run with r0 < 0.; Aborting");
	
	getInputFloat (&inp, "chain_stiff", &tmpf, 1);
	_chain_stiff = (number) tmpf;
	if (_chain_stiff < 0.) throw oxDNAException ("(DirkInteractionBias.cpp) Can't run with chain_stiff < 0.; Aborting");

	getInputFloat (&inp, "chain_r0", &tmpf, 1);
	_chain_r0 = (number) tmpf;
	if (_chain_r0 < 0.) throw oxDNAException ("(DirkInteractionBias.cpp) Can't run with chain_r0 < 0.; Aborting");
	
	OX_LOG(Logger::LOG_INFO, "(DirkInteractionBias.cpp) Initializing Dirk's interaction with length %g, DHS_radius = %g, DHS_eps = %g, DHS_cutoff = %g, overall rcut = %g, hard_rcut %g, bias K %g, bias r0 %g", _length, _DHS_radius, _DHS_eps, _DHS_rcut, this->_rcut, _hard_rcut, _K, _r0);
}

template<typename number>
void DirkInteractionBias<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
number DirkInteractionBias<number>::pair_interaction (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->n3 == q || p->n5 == q) return pair_interaction_bonded (p, q, r, update_forces);
	else return pair_interaction_nonbonded (p, q, r, update_forces);
}

template<typename number>
number DirkInteractionBias<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (p == P_VIRTUAL || q == P_VIRTUAL) return (number) 0.;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	//return (number) 0.;
	number energy = _hard (p, q, r, update_forces);
	energy += _dipolar (p, q, r, update_forces);
	energy += _chain (p, q, r, update_forces);
	
	return energy;
}

template<typename number>
number DirkInteractionBias<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (p->n3 == q || p->n5 == q) return (number) 0.;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	if(r->norm() >= this->_sqr_rcut) return (number) 0;

	number energy = _hard(p, q, r, update_forces);
	energy += _dipolar (p, q, r, update_forces); 

	return energy;
}

template<typename number>
number DirkInteractionBias<number>::_hard(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	number rnorm = (*r).norm(); 
	
	// overlap between DHS
	LR_vector<number> pdhs = p->orientation.v3 * _length / 2.;
	LR_vector<number> qdhs = q->orientation.v3 * _length / 2.;
	LR_vector<number> drdhs = *r + qdhs - pdhs;
	number drdhsnorm = drdhs.norm();

	if (rnorm < _hard_sqr_rcut) {
		if (drdhsnorm <= (4. * _DHS_radius * _DHS_radius)) {
			this->set_is_infinite(true);
			return (number) 1.e12;
		}
	
		// overlap between cylinders
		if (rnorm < _length * _length + 0.5 * 0.5 + 0.01) {
			bool cylinder_overlap = InteractionUtils::cylinder_overlap (p, q, *r, _length);
			if (cylinder_overlap) {
				this->set_is_infinite (true);
				return (number) 1.e12;
			}
		}
		
		// helpers
		number diag_cyl = sqrt((1. + _length * _length) / 4.);
		number compar = _DHS_radius * _DHS_radius + diag_cyl * diag_cyl + 2. * _DHS_radius * diag_cyl;
		
		// dhsq - cylinderp
		LR_vector<number> drcyl = *r + qdhs;
		if (drcyl.norm () < compar) {
			if (InteractionUtils::cylinder_sphere_overlap (drcyl, p->orientation.v3, _length, _DHS_radius)) {
				this->set_is_infinite(true);
				return (number) 1.e12;
			}
		}
	
		// dhsp - cylinderq
		drcyl = -(*r - pdhs); 
		if (drcyl.norm () < compar) {
			if (InteractionUtils::cylinder_sphere_overlap (drcyl, q->orientation.v3, _length, _DHS_radius)) {
				this->set_is_infinite(true);
				return (number) 1.e12;
			}
		}
	}
	return (number) 0.;
}

template<typename number>
number DirkInteractionBias<number>::_dipolar (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	
	LR_vector<number> pdhs = p->orientation.v3 * _length / 2.;
	LR_vector<number> qdhs = q->orientation.v3 * _length / 2.;
	LR_vector<number> drdhs = *r + qdhs - pdhs;
	number drdhsnorm = drdhs.norm();
	
	if (drdhsnorm > _DHS_sqr_rcut) return (number) 0.f;
	
	// direct part
	number energy = - (1. / (drdhsnorm * sqrt(drdhsnorm))) * (3.f * drdhs.z * drdhs.z / drdhsnorm - (number) 1.f);
	energy += - _DHS_rf_fact; 
	return energy;
}

template<typename number>
number DirkInteractionBias<number>::_chain (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	number energy = (number) 0.;
	if (p->n3 == q || p->n5 == q) {
		LR_vector<number> pdhs = p->orientation.v3 * _length / 2.;
		LR_vector<number> qdhs = q->orientation.v3 * _length / 2.;
		LR_vector<number> drdhs = *r + qdhs - pdhs;
		number dr = sqrt(drdhs.norm());
		energy += _chain_stiff * (dr *dr - 2. * dr * _chain_r0 + _chain_r0 * _chain_r0) / (number) 2.;
	}
	return energy;
}

template<typename number>
number DirkInteractionBias<number>::_bias (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	
	return (number) 0.f;
}

template<typename number>
void DirkInteractionBias<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
void DirkInteractionBias<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
bool DirkInteractionBias<number>::generate_random_configuration_overlap (BaseParticle<number> *p, BaseParticle<number> *q, number box_side) {
	LR_vector<number> dr = q->pos.minimum_image (p->pos, box_side);
	number energy = _dirk_pot (p, q, &dr, false);  
	this->set_is_infinite (false);
	return (energy > (number) 100.f);
}

template<typename number>
void DirkInteractionBias<number>::read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_strands, particles);
	int my_N, my_N_strands;

	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if(!topology.good()) throw oxDNAException("(DirkInteractionBias.cpp) Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	int strand, i = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#') continue;
		if(i == N_from_conf) throw oxDNAException("(DirkInteractionBias.cpp) Too many particles found in the topology file (should be %d), aborting", N_from_conf);

		int tmpn3, tmpn5;
		int res = sscanf(line, "%d %d %d", &strand, &tmpn3, &tmpn5);

		if(res < 3) throw oxDNAException("(DirkInteractionBias.cpp) Line %d of the topology file has an invalid syntax. Expectiong <int> <int> <int> (chain id, left neighbour id, right neighbour id)", i+2);

		BaseParticle<number> *p = particles[i];

		if(tmpn3 < 0) p->n3 = P_VIRTUAL;
		else p->n3 = particles[tmpn3];
		if(tmpn5 < 0) p->n5 = P_VIRTUAL;
		else p->n5 = particles[tmpn5];

		// store the strand id
		// for a design inconsistency, in the topology file
		// strand ids start from 1, not from 0
		p->strand_id = strand - 1;

		p->index = i;
		i++;
	}

	if(i < N_from_conf) throw oxDNAException ("(DirkInteractionBias.cpp) Not enough particles found in the topology file (found %d, should be %d). Aborting", i, N_from_conf);

	topology.close();

	if(my_N != N_from_conf) throw oxDNAException ("(DirkInteractionBias.cpp) Number of lines in the configuration file and\nnumber of particles in the topology files don't match. Aborting");

	*N_strands = my_N_strands;
}

template class DirkInteractionBias<float>;
template class DirkInteractionBias<double>;


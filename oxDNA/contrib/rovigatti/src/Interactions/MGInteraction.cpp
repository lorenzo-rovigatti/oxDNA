#include "MGInteraction.h"

#include "Particles/CustomParticle.h"
#include <fstream>

template<typename number>
MGInteraction<number>::MGInteraction() :
				BaseInteraction<number, MGInteraction<number> >() {
	this->_int_map[BONDED] = &MGInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &MGInteraction<number>::pair_interaction_nonbonded;

	_MG_n = 0;
	_MG_alpha = _MG_beta = _MG_gamma = 0.;
}

template<typename number>
MGInteraction<number>::~MGInteraction() {

}

template<typename number>
void MGInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "MG_rfene", &_rfene, 1);
	getInputNumber(&inp, "MG_rcut", &this->_rcut, 1);
	getInputNumber(&inp, "MG_alpha", &_MG_alpha, 0);
	getInputInt(&inp, "MG_n", &_MG_n, 1);
	getInputString(&inp, "MG_bond_file", _bond_filename, 1);
}

template<typename number>
void MGInteraction<number>::init() {
	_sqr_rfene = SQR(_rfene);

	number rep_rcut = pow(2., 1. / _MG_n);
	_MG_sqr_rep_rcut = SQR(rep_rcut);
	this->_sqr_rcut = SQR(this->_rcut);

	OX_LOG(Logger::LOG_INFO, "MG: repulsive rcut: %lf (%lf), total rcut: %lf (%lf)", rep_rcut, _MG_sqr_rep_rcut, this->_rcut, this->_sqr_rcut);

	if(_MG_alpha != 0) {
		if(_MG_alpha < 0.) throw oxDNAException("MG_alpha may not be negative");
		_MG_gamma = M_PI / (2.25 - pow(2., 1. / 3.));
		_MG_beta = 2 * M_PI - 2.25 * _MG_gamma;
		OX_LOG(Logger::LOG_INFO, "MG: alpha = %lf, beta = %lf, gamma = %lf", _MG_alpha, _MG_beta, _MG_gamma);
	}
}

template<typename number>
number MGInteraction<number>::_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();

	if(sqr_r > _sqr_rfene) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		this->set_is_infinite(true);
		return 10e10;
	}

	number energy = -15 * _sqr_rfene * log(1. - sqr_r / _sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30 * _sqr_rfene / (_sqr_rfene - sqr_r);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number MGInteraction<number>::_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();
	if(sqr_r > this->_sqr_rcut) return (number) 0.;

	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0;
	number energy = 0;
	// cut-off for all the repulsive interactions
	if(sqr_r < _MG_sqr_rep_rcut) {
		number part = 1.;
		for(int i = 0; i < _MG_n / 2; i++) part /= sqr_r;
		energy += 4. * (part * (part - 1.)) + 1. - _MG_alpha;
		if(update_forces) force_mod += 4. * _MG_n * part * (2. * part - 1.) / sqr_r;
	}
	// attraction
	else if(_MG_alpha != 0.) {
		energy += 0.5 * _MG_alpha * (cos(_MG_gamma * sqr_r + _MG_beta) - 1.0);
		if(update_forces) force_mod += _MG_alpha * _MG_gamma * sin(_MG_gamma * sqr_r + _MG_beta);
	}

	if(update_forces) {
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number MGInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number MGInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = (number) 0.f;

	if(p->is_bonded(q)) {
		LR_vector<number> computed_r;
		if(r == NULL) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = this->_box->min_image(p->pos, q->pos);
				r = &computed_r;
			}
		}

		energy = _fene(p, q, r, update_forces);
		energy += _nonbonded(p, q, r, update_forces);
	}

	return energy;
}

template<typename number>
number MGInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	LR_vector<number> computed_r;
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _nonbonded(p, q, r, update_forces);
}

template<typename number>
void MGInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
		for(typename std::set<CustomParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++) {

		}
	}
}

template<typename number>
void MGInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new CustomParticle<number>();
	}
}

template<typename number>
int MGInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);

	int N;
	sscanf(line, "%d\n", &N);

	return N;
}

template<typename number>
void MGInteraction<number>::read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_strands, particles);

	char line[512];
	std::ifstream bond_file(_bond_filename.c_str());

	if(!bond_file.good()) throw oxDNAException("Can't read bond file '%s'. Aborting", _bond_filename.c_str());
	// skip the headers and the particle positions
	int to_skip = N_from_conf + 2;
	for(int i = 0; i < to_skip; i++) {
		bond_file.getline(line, 512);
	}
	if(!bond_file.good()) throw oxDNAException("The bond file '%s' does not contain the right number of lines. Aborting", _bond_filename.c_str());

	for(int i = 0; i < N_from_conf; i++) {
		int idx, n_bonds;
		bond_file >> idx >> n_bonds;
		idx--;
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[idx]);
		p->type = P_A;
		p->btype = n_bonds;
		p->strand_id = 0;
		p->n3 = p->n5 = P_VIRTUAL;
		for(int j = 0; j < n_bonds; j++) {
			int n_idx;
			bond_file >> n_idx;
			n_idx--;
			CustomParticle<number> *q = static_cast<CustomParticle<number> *>(particles[n_idx]);
			p->add_bonded_neigh(q);
		}
		if(i != idx) throw oxDNAException("There is something wrong with the bond file. Expected index %d, found %d\n", i, idx);
	}

	*N_strands = N_from_conf;
}

extern "C" MGInteraction<float> *make_MGInteraction_float() {
	return new MGInteraction<float>();
}

extern "C" MGInteraction<double> *make_MGInteraction_double() {
	return new MGInteraction<double>();
}

template class MGInteraction<float> ;
template class MGInteraction<double> ;


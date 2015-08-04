#include <fstream>

#include "TSPInteraction.h"

template<typename number>
TSPInteraction<number>::TSPInteraction() : BaseInteraction<number, TSPInteraction<number> >() {
	this->_int_map[BONDED] = &TSPInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &TSPInteraction<number>::pair_interaction_nonbonded;

	_alpha = NULL;
	_N_arms = _N_monomer_per_arm = NULL;
	_TSP_n = _N_stars = 0;

	_attractive_anchor = false;
	_only_chains = false;
	_only_intra = false;
}

template<typename number>
TSPInteraction<number>::~TSPInteraction() {
	if(_alpha != NULL) delete[] _alpha;
	if(_N_arms != NULL) delete[] _N_arms;
	if(_N_monomer_per_arm != NULL) delete[] _N_monomer_per_arm;
}

template<typename number>
void TSPInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "TSP_rfene", &_rfene, 1);

	_rfene_anchor = _rfene * 3;
	getInputNumber(&inp, "TSP_rfene_anchor", &_rfene_anchor, 0);

	getInputNumber(&inp, "TSP_lambda", &_TSP_lambda, 1);
	if(_TSP_lambda < 0.) {
		OX_LOG(Logger::LOG_INFO, "Invalid TSP_lambda value detected, setting it to 0");
		_TSP_lambda = 0.;
	}

	getInputNumber(&inp, "TSP_rcut", &this->_rcut, 1);

	getInputInt(&inp, "TSP_n", &_TSP_n, 1);

	getInputBool(&inp, "TSP_attractive_anchor", &_attractive_anchor, 0);
	getInputBool(&inp, "TSP_only_chains", &_only_chains, 0);
	getInputBool(&inp, "TSP_only_intra", &_only_intra, 0);
}

template<typename number>
void TSPInteraction<number>::init() {
	_sqr_rfene = SQR(_rfene);
	_sqr_rfene_anchor = SQR(_rfene_anchor);

	number rep_rcut = pow(2., 1./_TSP_n);
	_TSP_sqr_rep_rcut = SQR(rep_rcut);
	OX_LOG(Logger::LOG_INFO, "TSP: repulsive rcut: %lf (%lf)", rep_rcut, _TSP_sqr_rep_rcut);
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
number TSPInteraction<number>::_fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	bool anchor = ((TSPParticle<number> *)p)->is_anchor() || ((TSPParticle<number> *)q)->is_anchor();

	number sqr_r = r->norm();
	number sqr_rfene = (anchor && !_only_chains) ? _sqr_rfene_anchor : _sqr_rfene;

	if(sqr_r > sqr_rfene) {
		if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r), sqrt(_sqr_rfene));
		else return 1e10;
	}

	number energy = -15 * sqr_rfene * log(1 - sqr_r/sqr_rfene);

	if(update_forces) {
		// this number is the module of the force over r, so we don't have to divide the distance
		// vector by its module
		number force_mod = -30 * sqr_rfene / (sqr_rfene - sqr_r);
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number TSPInteraction<number>::_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	int int_type = q->type + p->type;

	number sqr_r = r->norm();
	// cut-off for the telechelic monomers
	if(sqr_r > this->_sqr_rcut) return (number) 0.;
	if(sqr_r < 0.5 && update_forces) {
		if(update_forces) {
			number force_mod = 1000;
			p->force -= *r * force_mod;
			q->force += *r * force_mod;
		}
		
		return 10000000.;
	}

	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0;
	number energy = 0;
	// cut-off for all the repulsive interactions
	if(sqr_r < _TSP_sqr_rep_rcut) {
		number part = pow(1./sqr_r, _TSP_n/2.);
		energy += 4 * (part * (part - 1.)) + 1.;
		if(update_forces) force_mod += 4 * _TSP_n * part * (2*part - 1) / sqr_r;
	}

	// telechelic monomers
	if(int_type == 2) {
		if(sqr_r < _TSP_sqr_rep_rcut) energy -= _TSP_lambda;
		// same as before except for a lambda in front of both energy and force
		else {
			number part = pow(1./sqr_r, _TSP_n/2.);
			energy += 4 * _TSP_lambda * part * (part - 1.);

			if(update_forces) force_mod += 4 * _TSP_lambda * _TSP_n * part * (2*part - 1) / sqr_r;
		}
	}

	if(update_forces) {
		p->force -= *r * force_mod;
		q->force += *r * force_mod;
	}

	return energy;
}

template<typename number>
number TSPInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number TSPInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = (number) 0.f;

	// if q == P_VIRTUAL we have to compute the bonded interactions acting between p and all its bonded neighbours
	if(q == P_VIRTUAL) {
		TSPParticle<number> *TSPp = (TSPParticle<number> *) p;
		for(typename set<TSPParticle<number> *>::iterator it = TSPp->bonded_neighs.begin(); it != TSPp->bonded_neighs.end(); it++) {
			energy += pair_interaction_bonded(p, *it, r, update_forces);
		}
	}
	else if(p->is_bonded(q)){
		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}

		if(p->index < q->index && update_forces) return energy;
		energy = _fene(p, q, r, update_forces);
		energy += _nonbonded(p, q, r, update_forces);
	}
	return energy;
}

template<typename number>
number TSPInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	// if p and q don't belong to the same star AND the user enabled the TSP_only_intra option then
	// the interaction should be 0
	if(_only_intra && p->strand_id != q->strand_id) return 0.;

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	if(r->norm() >= this->_sqr_rcut) return (number) 0;

	return _nonbonded(p, q, r, update_forces);
}

template<typename number>
void TSPInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {	
		TSPParticle<number> *p = (TSPParticle<number> *)particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);

		// check that the distance between bonded neighbor doesn't exceed a reasonable threshold
		number mind = 0.5;
		number maxd = _rfene;
		if(p->n3 != P_VIRTUAL && p->n3->n5 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n3;
			LR_vector<number> rv = p->pos - q->pos;
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n3->index, r);
		}

		if(p->n5 != P_VIRTUAL && p->n5->n3 != P_VIRTUAL) {
			BaseParticle<number> *q = p->n5;
			LR_vector<number> rv = p->pos - q->pos;
			number r = sqrt(rv*rv);
			if(r > maxd || r < mind)
				throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf)", i, p->n5->index, r);
		}
	}
}

template<typename number>
void TSPInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new TSPParticle<number>();
}

template<typename number>
int TSPInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);

	int my_N_stars;
	sscanf(line, "%d\n", &my_N_stars);

	int N_arms, N_monomer_per_arm;
	// because we have at least a particle per star (the anchor)
	int N_from_topology = my_N_stars;
	for(int i = 0; i < my_N_stars; i++) {
		if(!topology.good()) throw oxDNAException("Not enough stars found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_stars);
		topology.getline(line, 512);
		if(sscanf(line, "%d %d %*f\n", &N_arms, &N_monomer_per_arm) != 2) throw oxDNAException("The topology file does not contain any info on star n.%d", i);

		N_from_topology += N_monomer_per_arm * N_arms;
	}

	return N_from_topology;
}

template<typename number>
void TSPInteraction<number>::read_topology(int N_from_conf, int *N_stars, BaseParticle<number> **particles) {
	IBaseInteraction<number>::read_topology(N_from_conf, N_stars, particles);
	int my_N_stars;
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);

	topology.getline(line, 512);
	sscanf(line, "%d\n", &my_N_stars);

	_N_arms = new int[my_N_stars];
	_N_monomer_per_arm = new int[my_N_stars];
	_alpha = new number[my_N_stars];
	// because we have at least _N_stars anchors
	int N_from_topology = my_N_stars;
	for(int i = 0; i < my_N_stars; i++) {
		if(!topology.good()) throw oxDNAException("Not enough stars found in the topology file. There are only %d lines, there should be %d, aborting", i, my_N_stars);
		topology.getline(line, 512);
		float tmp;
		sscanf(line, "%d %d %f\n", _N_arms + i, _N_monomer_per_arm + i, &tmp);
		_alpha[i] = tmp;

		N_from_topology += _N_monomer_per_arm[i] * _N_arms[i];
	}

	// construct the topology, i.e. assign the right FENE neighbours to all the particles
	int p_ind = 0;
	for(int ns = 0; ns < my_N_stars; ns++) {
		int attractive_from = (int) round(_N_monomer_per_arm[ns] * (1. - _alpha[ns]));
		OX_DEBUG("Adding a TSP with %d arms, %d monomers per arm (of which %d repulsive)", _N_arms[ns], _N_monomer_per_arm[ns], attractive_from);
		TSPParticle<number> *anchor = (TSPParticle<number> *) particles[p_ind];
		anchor->flag_as_anchor();
		// this is an anchor: it has n_arms FENE neighbours since it is attached to each arm
		anchor->btype = anchor->type = (_attractive_anchor) ? P_B : P_A;
		anchor->n3 = anchor->n5 = P_VIRTUAL;
		anchor->strand_id = ns;

		_anchors.push_back(anchor);

		p_ind++;
		for(int na = 0; na < _N_arms[ns]; na++) {
			for(int nm = 0; nm < _N_monomer_per_arm[ns]; nm++) {
				int type = (nm >= attractive_from) ? P_B : P_A;
				TSPParticle<number> *p = (TSPParticle<number> *) particles[p_ind];
				p->type = p->btype = type;
				p->strand_id = ns;
				p->n3 = p->n5 = P_VIRTUAL;
				p->set_arm(na);
				// if it's the first monomer of the arm then add also the anchor
				// as a bonded neighbour
				if(nm == 0) {
					p->add_bonded_neigh(anchor);
					p->n3 = anchor;
				}
				else {
					TSPParticle<number> *p_n = (TSPParticle<number> *) particles[p_ind-1];
					p->n3 = p_n;
					p_n->n5 = p;
					p->add_bonded_neigh(p_n);
				}

				p_ind++;
			}
		}
	}

	if(N_from_topology < N_from_conf) throw oxDNAException ("Not enough particles found in the topology file (should be %d). Aborting", N_from_conf);

	topology.close();

	*N_stars = my_N_stars;
	_N_stars = my_N_stars;
}

template<typename number>
bool TSPInteraction<number>::_insert_anchor(BaseParticle<number> **particles, BaseParticle<number> *p, Cells<number> *c) {
	int i = 0;
	bool inserted = false;

	do {
		p->pos = LR_vector<number>(drand48()*this->_box_side, drand48()*this->_box_side, drand48()*this->_box_side);
		inserted = !_does_overlap(particles, p, c);
		i++;
	} while(!inserted && i < MAX_INSERTION_TRIES);

	p->orientation = Utils::get_random_rotation_matrix_from_angle<number> (M_PI);

	return (i != MAX_INSERTION_TRIES);
}

template<typename number>
bool TSPInteraction<number>::_does_overlap(BaseParticle<number> **particles, BaseParticle<number> *p, Cells<number> *c) {
	// here we take into account the non-bonded interactions
	vector<BaseParticle<number> *> neighs = c->get_neigh_list(p, true);
	for(unsigned int n = 0; n < neighs.size(); n++) {
		BaseParticle<number> *q = neighs[n];
		// particles with an index larger than p->index have not been inserted yet
		if(q->index < p->index && this->generate_random_configuration_overlap(p, q, this->_box_side)) return true;
	}

	return false;
}

template<typename number>
void TSPInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N, number box_side) {
	Cells<number> c(N, this->_box);
	c.init(particles, this->_rcut);

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		p->pos = LR_vector<number>(0., 0., 0.);
	}

	c.global_update();

	int straight_till = 20;

	int i = 0;
	number max_dist = 0.;

	if(_N_stars > 1 && !_only_chains) OX_LOG(Logger::LOG_WARNING, "Can't reliably generate TSP configurations containing more than one star. Be careful.");

	for(int ns = 0; ns < _N_stars; ns++) {
		BaseParticle<number> *anchor = particles[i];
		anchor->index = i;
		if(!_insert_anchor(particles, anchor, &c)) throw oxDNAException("Can't insert particle number %d", i);
		c.single_update(anchor);
		i++;

		LR_vector<number> anchor_pos = anchor->pos;

		for(int na = 0; na < _N_arms[ns]; na++) {
			LR_vector<number> dir = Utils::get_random_vector<number>();
			number arm_dist = 1.;
			// if we simulate TSPs then the distance between an anchor and the first monomer of each arm
			// should be more than one sigma to avoid overlaps
			if(!_only_chains) arm_dist += 1.5;
			// if we simulate chains then we want monomers to be more than one sigma apart
			else dir *= 1.01;

			LR_vector<number> last_pos = anchor_pos;
			for(int nm = 0; nm < _N_monomer_per_arm[ns]; nm++) {
				BaseParticle<number> *p = particles[i];
				p->pos = last_pos;
				if(nm == 0) p->pos += dir*arm_dist;
				else {
					if(nm < straight_till && !_only_chains) p->pos += dir;
					else {
						number angle = drand48() * (2*M_PI-0.2) + 0.2;
						LR_matrix<number> R = Utils::get_random_rotation_matrix_from_angle<number>(angle);
						p->pos += R*dir;
					}
				}

				// if there are only chains then we want to generate non-straight chains
				// and hence we extract a new random direction
				if(_only_chains) dir = Utils::get_random_vector<number>();

				p->index = i;
				// if we simulate TSPs then we check for overlaps only for the first monomer
				if(nm == 0 && !_only_chains) {
					if(_does_overlap(particles, p, &c)) {
						na--;
						break;
					}
				}
				// if we simulate chains, on the other hand, we always check for overlaps
				if((_only_chains || nm >= straight_till) && _does_overlap(particles, p, &c)) {
					nm--;
					continue;
				}

				p->orientation = Utils::get_random_rotation_matrix_from_angle<number>(M_PI);
				c.single_update(p);
				last_pos = p->pos;
				number dist = sqrt(p->pos.sqr_distance(anchor_pos));
				if(dist > max_dist) max_dist = dist;

				i++;
			}
		}
	}
	OX_LOG(Logger::LOG_INFO, "Max distance: %f", max_dist);
}

template class TSPInteraction<float>;
template class TSPInteraction<double>;


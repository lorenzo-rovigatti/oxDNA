/**
 * @file    VMMC.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 * TODO
 * * fix acceptance
 * * handle diverging interaction
 * * record statistics
 */

#include "VMMC.h"
#include "../../Lists/Cells.h"

#define VMMC_ROTATION (1)
#define VMMC_TRANSLATION (2)

/// traslation
template<typename number>
VMMC<number>::VMMC (ConfigInfo<number> * Info) : BaseMove<number>(Info) {
	_max_move_size_sqr = -1.;
	_max_move_size = 1.;
	_max_cluster_size = -1;
}

template<typename number>
VMMC<number>::~VMMC () {
	delete _particles_old;
}

template<typename number>
void VMMC<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta_tras", &_delta_tras, 1);
	getInputNumber (&inp, "delta_rot", &_delta_rot, 1);
	getInputNumber (&inp, "prob", &this->prob, 0);
	
	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			throw oxDNAException ("Cannot run MC2 VMMC with verlet lists. Set list_type=cells in the input file");
		}
	}
	getInputInt(&inp, "max_cluster_size", &_max_cluster_size, 0);
	getInputNumber(&inp, "max_move_size", &_max_move_size, 0);

}

template<typename number>
void VMMC<number>::init() {
	BaseMove<number>::init();

	// setting the maximum displacement
	_max_move_size_sqr = _max_move_size * _max_move_size;
	
	// fix maxclust if evidently wrong
	if (_max_cluster_size < 1) {
		_max_cluster_size = *this->_Info->N;
	}
	if (_max_cluster_size > *this->_Info->N) {
		OX_LOG(Logger::LOG_WARNING, "(VMMC.cpp) maxclust > N does not make sense, setting it to N = %i", this->_Info->N);
		_max_cluster_size = *this->_Info->N;
	}

	// here we only use this array as a storage place for positions and orientations;
	// the topology is not set for _particles_old
	_particles_old = new BaseParticle<number>*[*this->_Info->N]; 
	for (int i = 0; i < *this->_Info->N; i ++) _particles_old[i] = new BaseParticle<number>();

	OX_LOG(Logger::LOG_INFO, "(VMMC.cpp) VMMC move initialized with T=%g, delta_tras=%g, delta_rot=%g, prob=%g, max_move_size=%g, max_cluster_size=%d", this->_T, _delta_tras, _delta_rot, this->prob, _max_move_size, _max_cluster_size);
	printf ("adjust: %d\n", this->_adjust_moves);

	if( ! dynamic_cast<Cells<number> *> (this->_Info->lists) )
			throw oxDNAException ("Cannot run MC2 VMMC with selected lists. Set list_type = cells in the input file");
	//throw oxDNAException ("Aborting at the end of VMMC::init for debugging purposes...");
	//abort();
}

// TODO: define store, restore, rerestore
template<typename number>
inline void VMMC<number>::_store_particle(BaseParticle<number> * src) {
	BaseParticle<number> * dst = _particles_old[src->index];

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->set_positions();
	dst->ext_potential = src->ext_potential;

	return;
}

template<typename number>
inline void VMMC<number>::_restore_particle (BaseParticle<number> * dst) {
	BaseParticle<number> *src = _particles_old[dst->index];

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->set_positions();
	dst->ext_potential = src->ext_potential;

	return;
}


template<typename number>
number VMMC<number>::build_cluster (movestr<number> * moveptr, int maxsize) {
	int nclust = 1;
	_clust.push_back(moveptr->seed);

	BaseParticle<number> * pp, * qq;

	set<int> prelinked_particles; //number of prelinked particles

	// CLUSTER GENERATION
	pp = this->_Info->particles[_clust[0]];
	pp->inclust = true;

	std::set<ParticlePair<number> > possible_links;

	typename vector<ParticlePair<number> >::iterator itb = pp->affected.begin();
	for(; itb != pp->affected.end(); itb++) {
		if (pp != itb->first && pp != itb->second) continue;
		if (pp == itb->first) qq = itb->second;
		else qq = itb->first;
		if (qq->inclust == false) possible_links.insert(ParticlePair<number>(pp, qq));
	}
	std::vector<BaseParticle<number> *> neighs1 = this->_Info->lists->get_neigh_list(pp);
	_store_particle(pp);
	_move_particle(moveptr, pp);
	this->_Info->lists->single_update(pp);
	std::vector<BaseParticle<number> *> neighs2 = this->_Info->lists->get_neigh_list(pp);
	_restore_particle (pp);
	this->_Info->lists->single_update(pp);
	typename std::vector<BaseParticle<number> *>::iterator itt = neighs1.begin();
	for (itt = neighs1.begin(); itt != neighs1.end(); itt ++) possible_links.insert(ParticlePair<number>(pp, *itt)); 
	for (itt = neighs2.begin(); itt != neighs2.end(); itt ++) possible_links.insert(ParticlePair<number>(pp, *itt)); 

	number  E_pp_moved, E_qq_moved, E_old;

	while (possible_links.size() > 0 && nclust <= maxsize) {
		// get a random link
		typename std::set<ParticlePair<number> >::iterator itc = possible_links.begin();

		itc = possible_links.begin();
		int n = rand() % ((int)possible_links.size());
		for (int kk = 0; kk < n; kk ++) itc ++;
		
		pp = (*itc).first;
		qq = (*itc).second;
		assert(pp->inclust || qq->inclust);

		// we remove this link and select a new one if both 
		// particles are in the cluster
		if (pp->inclust && qq->inclust) {
			possible_links.erase(ParticlePair<number>(pp, qq));
			continue;
		}

		// invertiamo se le particelle sono ordinate male
		// after this, now pp is in the cluster, qq is NOT in the cluster
		if (!pp->inclust) {
			BaseParticle<number> * tmp = qq;
			qq = pp;
			pp = tmp;
		}
		
		// now we check if pp and qq are bonded; if so, we store two auxiliary pointers
		// fpp and fqq to preserve the order of the particles when calling bonded interactions
		bool pp_qq_bonded = false;
		BaseParticle<number> * fpp = pp;
		BaseParticle<number> * fqq = qq;
		for(itb = pp->affected.begin(); itb != pp->affected.end(); itb++) {
			if ((*itb).first == pp && (*itb).second == qq) {
				fpp = pp;
				fqq = qq;
				pp_qq_bonded = true;		
			}
			if ((*itb).first == qq && (*itb).second == pp) {
				fpp = qq;
				fqq = pp;
				pp_qq_bonded = true;
			}
		}

		if (pp_qq_bonded) E_old = this->_Info->interaction->pair_interaction_bonded(fpp, fqq);
		else E_old = this->_Info->interaction->pair_interaction_nonbonded(pp, qq);
		assert (this->_Info->interaction->get_is_infinite() == false);

		_store_particle(pp);
		_move_particle(moveptr, pp);
		if (pp_qq_bonded) E_pp_moved = this->_Info->interaction->pair_interaction_bonded(fpp, fqq);
		else E_pp_moved = this->_Info->interaction->pair_interaction_nonbonded(pp, qq);
		_restore_particle(pp);
		bool force_prelink = this->_Info->interaction->get_is_infinite();
		this->_Info->interaction->set_is_infinite(false);

		number p1 = VMMC_link(E_pp_moved, E_old);
		if (p1 < (number)0.) p1 = (number) 0.f;

		if (force_prelink || p1 > this->_next_rand()) {
			_store_particle(qq);
			_move_particle(moveptr, qq);
			if (pp_qq_bonded) E_qq_moved = this->_Info->interaction->pair_interaction_bonded(fpp, fqq);
			else E_qq_moved = this->_Info->interaction->pair_interaction_nonbonded(pp, qq);
			_restore_particle(qq);
			force_prelink = this->_Info->interaction->get_is_infinite();
			this->_Info->interaction->set_is_infinite(false);

			number p2 = VMMC_link(E_qq_moved, E_old);
			if (p2 > (number)1.f) p2 = (number) 1.f;

			if (force_prelink || (p2 / p1) > this->_next_rand()) {
				// particle qq recruited
				qq->inclust = true;
				_clust.push_back(qq->index);
				nclust ++;
				
				for(itb = qq->affected.begin(); itb != qq->affected.end(); itb++) {
					BaseParticle<number> * tmp;
					if (qq != itb->first && qq != itb->second) continue;
					if (qq == itb->first) tmp = itb->second;
					else tmp = itb->first;
					if (tmp->inclust == false) possible_links.insert(ParticlePair<number>(qq, tmp));
				}

				neighs1 = this->_Info->lists->get_neigh_list(qq);
				_store_particle(qq);
				_move_particle(moveptr, qq);
				this->_Info->lists->single_update(qq);
				neighs2 = this->_Info->lists->get_neigh_list(qq);
				_restore_particle (qq);
				this->_Info->lists->single_update(qq);
				for (itt = neighs1.begin(); itt != neighs1.end(); itt ++) 
					if (!(*itt)->inclust)
						possible_links.insert(ParticlePair<number>(qq, *itt)); 
				for (itt = neighs2.begin(); itt != neighs2.end(); itt ++)
					if (!(*itt)->inclust)
						possible_links.insert(ParticlePair<number>(qq, *itt)); 
			}
			else {
				prelinked_particles.insert(qq->index);
			}
		}
		else {
			;
		}
		possible_links.erase(ParticlePair<number>(pp, qq));
	}

	// we break here if the cluster is too large
	if (nclust > maxsize) {
		return 0.;
	}

	for (int i = 0; i < nclust; i ++) {
		prelinked_particles.erase(_clust[i]);
	}
	if (prelinked_particles.size() > 0) {
		return (number) 0.;
	}

	for (int i = 0; i < nclust; i++) {
		pp = this->_Info->particles[_clust[i]];
		_store_particle(pp);
		_move_particle(moveptr, pp);
		if (pp->pos.sqr_distance(this->_particles_old[pp->index]->pos) > this->_max_move_size_sqr) {
			//printf("Cluster move done: exiting because the move is too large, setting pprime = 0.\n\n");
			return 0.;
		}
		this->_Info->lists->single_update(pp);
	}

	return 1.;
}


template<typename number>
void VMMC<number>::apply (llint curr_step) {
	this->_attempted += 1;

	// clear the cluster
	if (_clust.size() > 0) _clust.clear();

	// generate the move
	int pi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> *p = this->_Info->particles[pi];
	movestr<number> move;
	move.seed = pi;
	move.seed_strand_id = p->strand_id;
	//move.type = (drand48() < 0.5) ? VMMC_TRANSLATION : VMMC_ROTATION;
	move.type = VMMC_TRANSLATION;
	if (p->is_rigid_body() && (drand48() > 0.5)) move.type = VMMC_ROTATION;
	if (move.type == VMMC_TRANSLATION) {
		move.t = LR_vector<number> (Utils::gaussian<number>(), Utils::gaussian<number>(), Utils::gaussian<number>()) *_delta_tras;
	}
	else {
		// the translation vector is then interpreted as the point around which we rotate
		move.R = Utils::get_random_rotation_matrix_from_angle<number> (_delta_rot * Utils::gaussian<number>());
		move.Rt = (move.R).get_transpose();
		// WAS THIS move.t = this->_particles[move.seed]->int_centers[DNANucleotide<number>::BACK];
		// WE MAY WANT THIS move.t = this->_particles[move.seed]->int_centers[DNANucleotide<number>::BACK] + this->_particles[move.seed]->pos;
		move.t = this->_Info->particles[move.seed]->pos;
	}

	// build the cluster;
	number pprime;
	pprime = build_cluster (&move, _max_cluster_size);

	assert (_clust.size() >=1);
	int nclust = _clust.size();

	number delta_E_ext = 0.;

	// if we are not SURE to reject the move, we check the external
	// forces. Otherwise, there is no point.
	if (this->_Info->interaction->get_is_infinite() == false && pprime > 0.) {
		for (int l = 0; l < nclust; l++) {
			p = this->_Info->particles[_clust[l]];
			delta_E_ext += - p->ext_potential;
			p->set_ext_potential(curr_step, *this->_Info->box_side);
			delta_E_ext += + p->ext_potential;
		}
		pprime *= exp(-(1. / this->_T) * delta_E_ext);
	}

	if (this->_Info->interaction->get_is_infinite() == false && pprime > drand48()) {
		// move accepted
		this->_accepted += 1;
		
		// adjust move size
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			if (move.type == VMMC_TRANSLATION) {
				_delta_tras *= this->_acc_fact;
				if (_delta_tras > _max_move_size / 2.) _delta_tras = _max_move_size / 2.;
			}
			else {
				_delta_rot *= this->_acc_fact;
				if (_delta_rot > M_PI) _delta_rot = M_PI;
			}
		}
	}
	else {
		//move rejected
		for (int l = 0; l < nclust; l ++) {
			BaseParticle<number> * pp = this->_Info->particles[_clust[l]];
			_restore_particle (pp);
			this->_Info->lists->single_update(pp);
			pp->set_ext_potential(curr_step, *this->_Info->box_side);
		}
		this->_Info->interaction->set_is_infinite(false);

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			if (move.type == VMMC_TRANSLATION)_delta_tras /= this->_rej_fact;
			else _delta_rot /= this->_rej_fact;
		}
	}

	for (int l = 0; l < nclust; l ++) this->_Info->particles[_clust[l]]->inclust = false;

	return;
}

template<typename number>
inline void VMMC<number>::_move_particle(movestr<number> * moveptr, BaseParticle< number> *q) {
	if (moveptr->type == VMMC_TRANSLATION) {
		q->pos += moveptr->t;
	}
	else if (moveptr->type == VMMC_ROTATION) {
		LR_vector<number> dr;
		if (moveptr->seed_strand_id == q->strand_id) dr = q->pos - moveptr->t;
		else dr = q->pos.minimum_image(moveptr->t, *this->_Info->box_side);

		LR_vector<number> drp = moveptr->R * dr;
		q->pos = moveptr->t + drp;
		q->orientation = moveptr->R * q->orientation;
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();
	}
	else {
		;
	}

	return;
}


template class VMMC<float>;
template class VMMC<double>;

/*
 * MC_CPUBackend.cpp
 *
 *  Created on: 26/nov/2010
 *      Author: lorenzo
 */

#include "MC_CPUBackend.h"
#include "../Managers/SimManager.h"

template<typename number>
MC_CPUBackend<number>::MC_CPUBackend() : MCBackend<number>(), _particles_old(NULL) {
	this->_is_CUDA_sim = false;
	_enable_flip = false;
	_target_box = -1.;
	_box_tolerance = 1.e-8;
	_e_tolerance = -1.;
}

template<typename number>
MC_CPUBackend<number>::~MC_CPUBackend() {
	if (_particles_old != NULL) {
		// horrendous hack not to have forces deleted twice
		for(int i = 0; i < this->_N; i++) _particles_old[i]->N_ext_forces = 0;
		for(int i = 0; i < this->_N; i++) delete _particles_old[i];
		delete [] _particles_old;
	}
}

template<typename number>
void MC_CPUBackend<number>::get_settings(input_file &inp) {
	MCBackend<number>::get_settings(inp);

	getInputNumber<number> (&inp, "verlet_skin", &_verlet_skin, 0);
	getInputBool (&inp, "enable_flip", &_enable_flip, 0);
	if (_enable_flip) OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Enabling flip move");

	getInputNumber (&inp, "target_box", &_target_box, 0);
	number tmpf;
	if (getInputNumber (&inp, "target_volume", &tmpf, 0) == KEY_FOUND) _target_box = pow (tmpf, 1. / 3.);

	if (_target_box > 0.) getInputNumber (&inp, "box_tolerance", &_box_tolerance, 0);

	if (_target_box > 0.) {
		if (getInputNumber (&inp, "e_tolerance", &_e_tolerance, 0) == KEY_FOUND) {
			if (_e_tolerance < 0.) throw oxDNAException ("(MC_CPUBackend) Cannot run box adjustment with e_tolerance < 0.");
		}
	}
}

template<typename number>
void MC_CPUBackend<number>::init() {
	MCBackend<number>::init();

	_timer_move = TimingManager::instance()->new_timer(std::string("Rotations+Translations"), std::string("SimBackend"));
	if (this->_ensemble == MC_ENSEMBLE_NPT) _timer_box = TimingManager::instance()->new_timer(std::string("Volume Moves"), std::string("SimBackend"));
	_timer_lists = TimingManager::instance()->new_timer(std::string("Lists"));

	_particles_old = new BaseParticle<number>*[this->_N];
	this->_interaction->read_topology(this->_N, &this->_N_strands, _particles_old);
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = _particles_old[i];

		p->index = i;
		p->type = this->_particles[i]->type;
		p->init();

		for (int l = 0; l < this->_particles[i]->N_ext_forces; l ++) {
			p->add_ext_force(this->_particles[i]->ext_forces[l]);
		}
		p->copy_from(*this->_particles[i]);

		// this is needed for the first _compute_energy()
		this->_particles[i]->set_positions();
		this->_particles[i]->orientationT = this->_particles[i]->orientation.get_transpose();
	}

	// qui ci facciamo la lista delle interazioni bonded
	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		for (unsigned int n = 0; n < p->affected.size(); n ++) {
			BaseParticle<number> * p1 = p->affected[n].first;
			BaseParticle<number> * p2 = p->affected[n].second;
			number e = this->_interaction->pair_interaction_bonded (p1, p2);
			if (_stored_bonded_interactions.count(ParticlePair<number>(p1, p2)) == 0) _stored_bonded_interactions[ParticlePair<number>(p1,p2)] = e;
		}
	}

	// check that target_box makes sense and output details
	if (_target_box > 0.f) {
		if (_target_box < 2. * this->_interaction->get_rcut())
			throw oxDNAException("Cannot run box adjustment with target_box (%g) <  2 * rcut (2 * %g)", _target_box, this->_interaction->get_rcut()); 
		if (_e_tolerance < 0.f) _e_tolerance = 0.3 * this->_T * this->_N;	
		else _e_tolerance *= this->_N;

		OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Working to achieve taget_box=%g (Volume: %g), tolerance=%g, energy tolerance=%g (e/N)=%g", _target_box, pow(_target_box,3.), _box_tolerance, _e_tolerance, _e_tolerance / this->_N);
	}

	_compute_energy();
	if(this->_overlap == true) throw oxDNAException("There is an overlap in the initial configuration");
}

template<typename number>
inline number MC_CPUBackend<number>::_particle_energy(BaseParticle<number> *p, bool reuse) {
	number res = (number) 0.f;

	if(reuse) {
		// slightly better than a direct loop, since in this way we don't have to build
		// ParticlePair objects every time
		typename vector<ParticlePair<number> >::iterator it = p->affected.begin();
		for(; it != p->affected.end(); it++) res += _stored_bonded_interactions[*it];
	}
	else {
		typename vector<ParticlePair<number> >::iterator it = p->affected.begin();
		for(; it != p->affected.end(); it++) {
			number de = this->_interaction->pair_interaction_bonded(it->first, it->second);
			res += de;
			_stored_bonded_tmp[*it] = de;
		}
	}

	if(this->_interaction->get_is_infinite() == true) {
		this->_overlap = true;
		return (number) 1.e12;
	}

	std::vector<BaseParticle<number> *> neighs = this->_lists->get_neigh_list(p);
	for(unsigned int n = 0; n < neighs.size(); n++) {
		BaseParticle<number> *q = neighs[n];
		res += this->_interaction->pair_interaction_nonbonded(p, q);
		if(this->_interaction->get_is_infinite() == true) {
			this->_overlap = true;
			return (number) 1.e12;
		}
	}

	return res;
}

template<typename number>
void MC_CPUBackend<number>::_compute_energy() {
	this->_U = (number) 0;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		this->_U += this->_particle_energy(p);
	}

	this->_U *= (number) 0.5;
}

template<typename number>
inline void MC_CPUBackend<number>::_translate_particle(BaseParticle<number> *p) {
	p->pos.x += (drand48() - (number)0.5f)*this->_delta[MC_MOVE_TRANSLATION];
	p->pos.y += (drand48() - (number)0.5f)*this->_delta[MC_MOVE_TRANSLATION];
	p->pos.z += (drand48() - (number)0.5f)*this->_delta[MC_MOVE_TRANSLATION];
}

template<typename number>
inline void MC_CPUBackend<number>::_rotate_particle(BaseParticle<number> *p) {
	number t;
	LR_vector<number> axis;
	if (_enable_flip && drand48() < (1.f / 20.f)) {
		// flip move
		//fprintf (stderr, "Flipping %d\n", p->index);
		t = M_PI / 2.;
		axis = p->orientation.v1;
		axis.normalize();
	}
	else {
		// normal random move
		t = (drand48() - (number)0.5f) * this->_delta[MC_MOVE_ROTATION];
		axis = Utils::get_random_vector<number>();
	}

	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number)1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	p->orientation = p->orientation * R;
}

template<typename number>
void MC_CPUBackend<number>::sim_step(llint curr_step) {
	this->_mytimer->resume();

	for(int i = 0; i < this->_N; i++) {
		if (i > 0 && this->_interaction->get_is_infinite() == true) throw oxDNAException ("should not happen %d", i);
		if (this->_ensemble == MC_ENSEMBLE_NPT && drand48() < 1. / this->_N) {
			_timer_box->resume();
			// do npt move

			// useful for checking
			number oldE2 = this->_interaction->get_system_energy(this->_particles, this->_N, this->_lists);
			if (fabs((this->_U - oldE2)/oldE2) > 1.e-6) throw oxDNAException ("happened %g %g", this->_U, oldE2);
			if (fabs((this->_U - oldE2)/oldE2) > 1.e-6) printf ("### happened %g %g (%g) (%g)\n", this->_U, oldE2, 
					fabs(this->_U- oldE2), fabs((this->_U - oldE2)/oldE2));
			if (this->_interaction->get_is_infinite()) throw oxDNAException ("non ci siamo affatto");
			number oldE = this->_U;

			number old_box_side = this->_box_side;

			number dL = this->_delta[MC_MOVE_VOLUME] * (drand48() - (number) 0.5);

			this->_box_side += dL;
			this->_interaction->set_box_side(this->_box_side);
			this->_box->init(this->_box_side, this->_box_side, this->_box_side);
			this->_lists->change_box();

			number dExt = (number) 0.;
			for (int k = 0; k < this->_N; k ++) {
				BaseParticle<number> *p = this->_particles[k];
				dExt = -p->ext_potential;
				_particles_old[k]->pos = p->pos; 
				p->pos *= this->_box_side / old_box_side;
				p->set_ext_potential(curr_step, this->_box_side);
				dExt += -p->ext_potential;
			}
			//for (int i = 0; i < this->_N; i ++) this->_lists->single_update(this->_particles[i]);
			if(!this->_lists->is_updated()) {
				_timer_lists->resume();
				this->_lists->global_update();
				this->_N_updates++;
				_timer_lists->pause();
			}

			number newE = this->_interaction->get_system_energy(this->_particles, this->_N, this->_lists);
			number dE = newE - oldE + dExt;
			number V = this->_box_side * this->_box_side * this->_box_side;
			number Vold = old_box_side * old_box_side * old_box_side;
			number dV = V - Vold;

			this->_tries[MC_MOVE_VOLUME]++;

			bool second_factor;
			if (_target_box > (number) 0.f) {
				number V_target = _target_box * _target_box * _target_box;
				second_factor = fabs(V - V_target) < fabs (Vold - V_target) && dE < _e_tolerance;
			}
			else second_factor = exp (-(dE + this->_P * dV - this->_N * this->_T * log (V / Vold)) / this->_T) > drand48();

			if (this->_interaction->get_is_infinite() == false && second_factor) {
				// volume move accepted
				this->_accepted[MC_MOVE_VOLUME]++;
				this->_U = newE;
				if ((curr_step < this->_MC_equilibration_steps && this->_adjust_moves) || _target_box > 0.f) this->_delta[MC_MOVE_VOLUME] *= 1.03;

				_stored_bonded_interactions.clear();
				for (int k = 0; k < this->_N; k ++) {
					BaseParticle<number> *p = this->_particles[k];
					typename vector<ParticlePair<number> >::iterator it = p->affected.begin();
					for(; it != p->affected.end(); it++) {
						number e = this->_interaction->pair_interaction_bonded(it->first, it->second);
						if(_stored_bonded_interactions.count(*it) == 0) _stored_bonded_interactions[*it] = e;
					}
				}

			}
			else {
				// volume move rejected
				for (int k = 0; k < this->_N; k ++) {
					BaseParticle<number> *p = this->_particles[k];
					//p->pos /= this->_box_side / old_box_side;
					p->pos = _particles_old[k]->pos; 
					p->set_ext_potential(curr_step, this->_box_side);
				}
				this->_box_side = old_box_side;
				this->_interaction->set_box_side(this->_box_side);
				this->_box->init(this->_box_side, this->_box_side, this->_box_side);
				this->_lists->change_box();
				this->_interaction->set_is_infinite (false);
				//for (int i = 0; i < this->_N; i ++) this->_lists->single_update(this->_particles[i]);
				if(!this->_lists->is_updated()) {
					_timer_lists->resume();
					//printf ("updating lists after  rejection\n");
					this->_lists->global_update();
					this->_N_updates++;
					_timer_lists->pause();
				}
				if ((curr_step < this->_MC_equilibration_steps && this->_adjust_moves) || _target_box > 0.f) this->_delta[MC_MOVE_VOLUME] /= 1.01;
			}
			_timer_box->pause();
		}
		else {
			_timer_move->resume();
			// do normal move
			int pi = (int) (drand48() * this->_N);
			BaseParticle<number> *p = this->_particles[pi];

			int move = (drand48() < (number) 0.5f) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;
			if(!p->is_rigid_body()) move = MC_MOVE_TRANSLATION;

			this->_tries[move]++;
			number delta_E = -_particle_energy(p, true);
			//number delta_E = -_particle_energy(p, false);
			p->set_ext_potential (curr_step, this->_box_side);
			number delta_E_ext = -p->ext_potential;

			if(move == MC_MOVE_TRANSLATION) {
				_particles_old[pi]->pos = p->pos;
				_translate_particle(p);
			}
			else {
				_particles_old[pi]->orientation = p->orientation;
				_particles_old[pi]->orientationT = p->orientationT;
				_rotate_particle(p);
				p->orientationT = p->orientation.get_transpose();
				p->set_positions();
			}
			this->_lists->single_update(p);

			//get_time(&this->_timer, 2);
			if(!this->_lists->is_updated()) {
				//printf ("updating lists because of translations\n");
				_timer_lists->resume();
				this->_lists->global_update();
				this->_N_updates++;
				_timer_lists->pause();
			}
			//get_time(&this->_timer, 3);

			_stored_bonded_tmp.clear();
			delta_E += _particle_energy(p, false);
			p->set_ext_potential(curr_step, this->_box_side);
			delta_E_ext += p->ext_potential;

			// uncomment to check the energy at a given time step.
			// may be useful for debugging purposes
			//if (curr_step > 410000 && curr_step <= 420001)
			// printf("delta_E: %lf\n", (double)delta_E);

			if(!this->_overlap && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48())) {
				this->_accepted[move]++;
				this->_U += delta_E;
				if (curr_step < this->_MC_equilibration_steps && this->_adjust_moves) {
					this->_delta[move] *= 1.03;
					if (move == MC_MOVE_TRANSLATION && this->_delta[move] > _verlet_skin * sqrt(3.) / 2.)
						this->_delta[move] = _verlet_skin * sqrt(3.) / 2. ;
					if (move == MC_MOVE_ROTATION && this->_delta[move] > M_PI / 2.)
						this->_delta[move] = M_PI / 2.;
				}

				// slightly faster than doing a loop over the indexes of affected
				typename map<ParticlePair<number>, number>::iterator it = _stored_bonded_tmp.begin();
				for(; it != _stored_bonded_tmp.end(); it++) {
					_stored_bonded_interactions[it->first] = it->second;
				}

			}
			else {
				if(move == MC_MOVE_TRANSLATION) {
					p->pos = _particles_old[pi]->pos;
				}
				else {
					p->orientation = _particles_old[pi]->orientation;
					p->orientationT = _particles_old[pi]->orientationT;
					p->set_positions();
				}
				this->_lists->single_update(p);
				this->_interaction->set_is_infinite(false);
				if(curr_step < this->_MC_equilibration_steps && this->_adjust_moves) this->_delta[move] /= 1.01;

			}
			this->_overlap = false;
			_timer_move->pause();
		}
	}

	if (_target_box > 0.) {
		if (fabs(_target_box / this->_box_side - 1.) < _box_tolerance) {
			SimManager::stop = true;
			OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Box adjusted to %g (target %g, relative error %g)", this->_box_side, _target_box, fabs(_target_box / this->_box_side - 1.));
		}
	}

	this->_mytimer->pause();
}

template class MC_CPUBackend<float>;
template class MC_CPUBackend<double>;


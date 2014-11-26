/*
 * MC_CPUBackend.cpp
 *
 *  Created on: 26/nov/2010
 *      Author: lorenzo
 */

#include "MC_CPUBackend.h"

template<typename number>
MC_CPUBackend<number>::MC_CPUBackend() : MCBackend<number>(), _particles_old(NULL) {
	this->_is_CUDA_sim = false;
	// initialize the messages for the timings output
	this->_timer_msgs_number = 2;
	strncpy(this->_timer_msgs[0], "MC step", 256);
	strncpy(this->_timer_msgs[1], "Lists update", 256);
	_enable_flip = false;
}

template<typename number>
MC_CPUBackend<number>::~MC_CPUBackend() {
	if (_particles_old != NULL) {
		// horrendous hack not to have forces deleted twice
		for(int i = 0; i < this->_N; i++) _particles_old[i]->N_ext_forces = 0;
		for(int i = 0; i < this->_N; i++) delete _particles_old[i];
		delete [] _particles_old;
	}
	if(this->_N_updates > 0) divide_given_timing(&this->_timer, 1, this->_N / (double) this->_N_updates);
}

template<typename number>
void MC_CPUBackend<number>::get_settings(input_file &inp) {
	MCBackend<number>::get_settings(inp);

	//float skin = 1.0f;
	//getInputFloat(&inp, "verlet_skin", &skin, 0);
	//_verlet_skin = skin;
	
	//int tmpi;
	//if (getInputBoolAsInt(&inp, "enable_flip", &tmpi, 0) == KEY_FOUND) {
	//	if (tmpi > 0) {
	//		OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Enabling flip move...");
	//		_enable_flip = true;
	//	}
	//}

	getInputNumber<number> (&inp, "verlet_skin", &_verlet_skin, 0);
	getInputBool (&inp, "enable_flip", &_enable_flip, 0);
	if (_enable_flip) OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Enabling flip move");
}

template<typename number>
void MC_CPUBackend<number>::init(char conf_filename[256]) {
	MCBackend<number>::init(conf_filename);

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

	_compute_energy();
	if(this->_overlap == true) throw oxDNAException("There is an overlap in the initial configuration");
}

template<typename number>
inline number MC_CPUBackend<number>::_particle_energy(BaseParticle<number> *p) {
	number res = (number) 0.f;
	if(p->n3 != P_VIRTUAL) res += this->_interaction->pair_interaction_bonded(p, p->n3);
	if(p->n5 != P_VIRTUAL) res += this->_interaction->pair_interaction_bonded(p, p->n5);

	if (this->_interaction->get_is_infinite() == true) {
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
	this->_U = this->_U_hydr = this->_U_stack = (number) 0;
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		if(p->n3 != P_VIRTUAL) this->_U += this->_interaction->pair_interaction_bonded(p, p->n3);
		if(p->n5 != P_VIRTUAL) this->_U += this->_interaction->pair_interaction_bonded(p, p->n5);

		std::vector<BaseParticle<number> *> neighs = this->_lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle<number> *q = neighs[n];
			this->_U += this->_interaction->pair_interaction_nonbonded(p, q);
		}
	}

	this->_U *= (number) 0.5;
	this->_U_hydr *= (number) 0.5;
	this->_U_stack *= (number) 0.5;
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
	LR_vector<number> tmp;

	get_time(&this->_timer, 0);

	for(int i = 0; i < this->_N; i++) {
		if (i > 0 && this->_interaction->get_is_infinite() == true) throw oxDNAException ("should not happen %d", i);
		if (this->_ensemble == MC_ENSEMBLE_NPT && drand48() < 1. / this->_N) {
			// do npt move

			// might be useful for checking
			//number oldE2 = this->_interaction->get_system_energy(this->_particles, this->_N, this->_box_side);
			//if (fabs((this->_U - oldE2)/oldE2) > 1.e-6) throw oxDNAException ("happened %g %g", this->_U, oldE2);

			number oldE = this->_U;

			number old_box_side = this->_box_side;

			number dL = this->_delta[MC_MOVE_VOLUME] * (drand48() - (number) 0.5);

			this->_box_side += dL;
			this->_interaction->set_box_side(this->_box_side);

			number dExt = (number) 0.;
			for (int k = 0; k < this->_N; k ++) {
				BaseParticle<number> *p = this->_particles[k];
				dExt = -p->ext_potential;
				p->pos *= this->_box_side / old_box_side;
				p->set_ext_potential(curr_step, this->_box_side);
				dExt += -p->ext_potential;
			}
			//for (int i = 0; i < this->_N; i ++) this->_lists->single_update(this->_particles[i]);
			if(!this->_lists->is_updated()) {
				this->_lists->global_update();
				this->_N_updates++;
			}

			number newE = this->_interaction->get_system_energy(this->_particles, this->_N, this->_box_side);
			number dE = newE - oldE + dExt;
			number V = this->_box_side * this->_box_side * this->_box_side;
			number Vold = old_box_side * old_box_side * old_box_side;
			number dV = V - Vold;

			this->_tries[MC_MOVE_VOLUME]++;

			if (this->_interaction->get_is_infinite() == false && exp (-(dE + this->_P * dV - this->_N * this->_T * log (V / Vold)) / this->_T) > drand48()) {
				// volume move accepted
				this->_accepted[MC_MOVE_VOLUME]++;
				this->_U = newE;
				if (curr_step < this->_MC_equilibration_steps && this->_adjust_moves) this->_delta[MC_MOVE_VOLUME] *= 1.03;
			}
			else {
				// volume move rejected
				for (int k = 0; k < this->_N; k ++) {
					BaseParticle<number> *p = this->_particles[k];
					p->pos /= this->_box_side / old_box_side;
					p->set_ext_potential(curr_step, this->_box_side);
				}
				this->_box_side = old_box_side;
				this->_interaction->set_box_side(this->_box_side);
				this->_interaction->set_is_infinite (false);
				//for (int i = 0; i < this->_N; i ++) this->_lists->single_update(this->_particles[i]);
				if(!this->_lists->is_updated()) {
					this->_lists->global_update();
					this->_N_updates++;
				}
				if (curr_step < this->_MC_equilibration_steps && this->_adjust_moves) this->_delta[MC_MOVE_VOLUME] /= 1.01;
			}
		}
		else {
			// do normal move
			int pi = (int) (drand48() * this->_N);
			BaseParticle<number> *p = this->_particles[pi];

			int move = (drand48() < (number) 0.5f) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;

			this->_tries[move]++;
			number delta_E = -_particle_energy(p);
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

			get_time(&this->_timer, 2);
			if(!this->_lists->is_updated()) {
				this->_lists->global_update();
				this->_N_updates++;
			}
			get_time(&this->_timer, 3);

			delta_E += _particle_energy(p);
			p->set_ext_potential(curr_step, this->_box_side);
			delta_E_ext += p->ext_potential;

			// uncomment to check the energy at a given time step.
			// may be useful for debugging purposes
			//if (curr_step > 410000 && curr_step <= 420001)
			// printf("delta_E: %lf\n", (double)delta_E);

			if(!this->_overlap && ((delta_E + delta_E_ext) < 0 ||
					       exp(-(delta_E + delta_E_ext) / this->_T) > drand48())) {
				this->_accepted[move]++;
				this->_U += delta_E;
				if (curr_step < this->_MC_equilibration_steps && this->_adjust_moves) {
					this->_delta[move] *= 1.03;
					if (move == MC_MOVE_TRANSLATION && this->_delta[move] > _verlet_skin * sqrt(3.) / 2.)
						this->_delta[move] = _verlet_skin * sqrt(3.) / 2. ;
					if (move == MC_MOVE_ROTATION && this->_delta[move] > M_PI / 2.)
						this->_delta[move] = M_PI / 2.;
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
				if (curr_step < this->_MC_equilibration_steps && this->_adjust_moves) this->_delta[move] /= 1.01;

			}
			this->_overlap = false;
		}
	}

	get_time(&this->_timer, 1);
	process_times(&this->_timer);
}

template class MC_CPUBackend<float>;
template class MC_CPUBackend<double>;


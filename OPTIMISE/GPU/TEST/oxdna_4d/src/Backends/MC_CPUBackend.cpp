/*
 * MC_CPUBackend.cpp
 *
 *  Created on: 26/nov/2010
 *      Author: lorenzo
 */

#include "MC_CPUBackend.h"

#include "../Particles/BaseParticle.h"
#include "../Observables/ObservableOutput.h"
#include "../Managers/SimManager.h"

MC_CPUBackend::MC_CPUBackend() :
				MCBackend() {
	_enable_flip = false;
	_target_box = -1.;
	_box_tolerance = 1.e-8;
	_e_tolerance = -1.;
	_verlet_skin = -1.;
}

MC_CPUBackend::~MC_CPUBackend() {
	for(auto p: _particles_old) {
		delete p;
	}
}

void MC_CPUBackend::get_settings(input_file &inp) {
	MCBackend::get_settings(inp);

	getInputNumber(&inp, "verlet_skin", &_verlet_skin, 0);
	getInputBool(&inp, "enable_flip", &_enable_flip, 0);
	if(_enable_flip) OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Enabling flip move");

	getInputNumber(&inp, "target_box", &_target_box, 0);
	number tmpf;
	if(getInputNumber(&inp, "target_volume", &tmpf, 0) == KEY_FOUND) _target_box = pow(tmpf, 1. / 3.);

	if(_target_box > 0.) getInputNumber(&inp, "box_tolerance", &_box_tolerance, 0);

	if(_target_box > 0.) {
		if(getInputNumber(&inp, "e_tolerance", &_e_tolerance, 0) == KEY_FOUND) {
			if(_e_tolerance < 0.) throw oxDNAException("(MC_CPUBackend) Cannot run box adjustment with e_tolerance < 0.");
		}
	}
}

void MC_CPUBackend::init() {
	MCBackend::init();

	_timer_move = TimingManager::instance()->new_timer(std::string("Rotations+Translations"), std::string("SimBackend"));
	if(_ensemble == MC_ENSEMBLE_NPT) _timer_box = TimingManager::instance()->new_timer(std::string("Volume Moves"), std::string("SimBackend"));
	_timer_lists = TimingManager::instance()->new_timer(std::string("Lists"));

	_particles_old.resize(N());
	_interaction->read_topology(&_N_strands, _particles_old);
	for(int i = 0; i < N(); i++) {
		BaseParticle *p = _particles_old[i];

		p->index = i;
		p->type = _particles[i]->type;
		p->init();

		for(auto ext_force : _particles[i]->ext_forces) {
			p->add_ext_force(ext_force);
		}
		p->copy_from(*_particles[i]);

		// this is needed for the first _compute_energy()
		_particles[i]->set_positions();
		_particles[i]->orientationT = _particles[i]->orientation.get_transpose();
	}

	// here we build the list of bonded interactions
	for(auto p: _particles) {
		for(unsigned int n = 0; n < p->affected.size(); n++) {
			BaseParticle * p1 = p->affected[n].first;
			BaseParticle * p2 = p->affected[n].second;
			number e = _interaction->pair_interaction_bonded(p1, p2);
			if(_stored_bonded_interactions.count(ParticlePair(p1, p2)) == 0) _stored_bonded_interactions[ParticlePair(p1, p2)] = e;
		}
	}

	// check that target_box makes sense and output details
	if(_target_box > 0.f) {
		if(_target_box < 2. * _interaction->get_rcut()) {
			throw oxDNAException("Cannot run box adjustment with target_box (%g) <  2 * rcut (2 * %g)", _target_box, _interaction->get_rcut());
		}
		if(_e_tolerance < 0.f) {
			_e_tolerance = 0.3 * N();
		}
		else {
			_e_tolerance *= N();
		}

		OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Working to achieve taget_box=%g (Volume: %g), tolerance=%g, energy tolerance=%g (e/N)=%g", _target_box, pow(_target_box,3.), _box_tolerance, _e_tolerance, _e_tolerance / N());
	}

	_compute_energy();
	if(_overlap == true) throw oxDNAException("(MC_CPUBackend) There is an overlap in the initial configuration");
}

inline number MC_CPUBackend::_particle_energy(BaseParticle *p, bool reuse) {
	number res = (number) 0.f;

	if(reuse) {
		// slightly better than a direct loop, since in this way we don't have to build
		// ParticlePair objects every time
		for(auto &pair : p->affected) {
			res += _stored_bonded_interactions[pair];
		}
	}
	else {
		for(auto &pair : p->affected) {
			number de = _interaction->pair_interaction_bonded(pair.first, pair.second);
			res += de;
			_stored_bonded_tmp[pair] = de;
		}
	}

	if(_interaction->get_is_infinite() == true) {
		_overlap = true;
		return (number) 1.e12;
	}

	std::vector<BaseParticle *> neighs = _lists->get_neigh_list(p);
	for(unsigned int n = 0; n < neighs.size(); n++) {
		BaseParticle *q = neighs[n];
		res += _interaction->pair_interaction_nonbonded(p, q);
		if(_interaction->get_is_infinite() == true) {
			_overlap = true;
			return (number) 1.e12;
		}
	}

	return res;
}

void MC_CPUBackend::_compute_energy() {
	_interaction->begin_energy_computation();

	_U = (number) 0;
	for(auto p: _particles) {
		_U += _particle_energy(p);
	}

	_U *= (number) 0.5;
}

inline void MC_CPUBackend::_translate_particle(BaseParticle *p) {
	p->pos.x += (drand48() - (number) 0.5f) * _delta[MC_MOVE_TRANSLATION];
	p->pos.y += (drand48() - (number) 0.5f) * _delta[MC_MOVE_TRANSLATION];
	p->pos.z += (drand48() - (number) 0.5f) * _delta[MC_MOVE_TRANSLATION];
}

inline void MC_CPUBackend::_rotate_particle(BaseParticle *p) {
	number t;
	LR_vector axis;
	if(_enable_flip && drand48() < (1.f / 20.f)) {
		// flip move
		//fprintf (stderr, "Flipping %d\n", p->index);
		t = M_PI / 2.;
		axis = p->orientation.v1;
		axis.normalize();
	}
	else {
		// normal random move
		t = (drand48() - (number) 0.5f) * _delta[MC_MOVE_ROTATION];
		axis = Utils::get_random_vector();
	}

	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number) 1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	p->orientation = p->orientation * R;
}

void MC_CPUBackend::sim_step() {
	_mytimer->resume();

	for(int i = 0; i < N(); i++) {
		if(i > 0 && _interaction->get_is_infinite() == true) {
			throw oxDNAException("should not happen %d", i);
		}
		if(_ensemble == MC_ENSEMBLE_NPT && drand48() < 1. / N()) {
			_timer_box->resume();
			// do npt move

			// useful for checking
			//number oldE2 = _interaction->get_system_energy(_particles, _N, _lists);
			//if (fabs((_U - oldE2)/oldE2) > 1.e-6) throw oxDNAException ("happened %g %g", _U, oldE2);
			//if (fabs((_U - oldE2)/oldE2) > 1.e-6) printf ("### happened %g %g (%g) (%g)\n", _U, oldE2,
			//		fabs(_U- oldE2), fabs((_U - oldE2)/oldE2));
			if(_interaction->get_is_infinite()) throw oxDNAException("non ci siamo affatto");
			number oldE = _U;
			number oldV = _box->V();

			LR_vector box_sides = _box->box_sides();
			LR_vector old_box_sides = box_sides;

			number dL = _delta[MC_MOVE_VOLUME] * (drand48() - (number) 0.5);

			// isotropic move
			box_sides.x += dL;
			box_sides.y += dL;
			box_sides.z += dL;

			_box->init(box_sides.x, box_sides.y, box_sides.z);

			number dExt = (number) 0.;
			for(int k = 0; k < N(); k++) {
				BaseParticle *p = _particles[k];
				dExt = -p->ext_potential;
				_particles_old[k]->pos = p->pos;
				p->pos.x *= box_sides[0] / old_box_sides[0];
				p->pos.y *= box_sides[1] / old_box_sides[1];
				p->pos.z *= box_sides[2] / old_box_sides[2];
				p->set_ext_potential(current_step(), _box.get());
				dExt += -p->ext_potential;
			}
			if(!_lists->is_updated()) {
				_timer_lists->resume();
				_lists->global_update();
				_N_updates++;
				_timer_lists->pause();
			}

			number newE = _interaction->get_system_energy(_particles, _lists.get());
			number dE = newE - oldE + dExt;
			number V = _box->V();
			number dV = V - oldV;

			_tries[MC_MOVE_VOLUME]++;

			bool second_factor;
			if(_target_box > (number) 0.f) {
				number V_target = _target_box * _target_box * _target_box;
				second_factor = fabs(V - V_target) < fabs(oldV - V_target) && dE < _e_tolerance;
			}
			else {
				second_factor = exp(-(dE + _P * dV - N() * _T * log(V / oldV)) / _T) > drand48();
			}

			if(_interaction->get_is_infinite() == false && second_factor) {
				// volume move accepted
				_accepted[MC_MOVE_VOLUME]++;
				_U = newE;
				if((current_step() < _MC_equilibration_steps && _adjust_moves) || _target_box > 0.f) _delta[MC_MOVE_VOLUME] *= 1.03;

				_stored_bonded_interactions.clear();
				for(auto p: _particles) {
					for(auto &pair : p->affected) {
						number e = _interaction->pair_interaction_bonded(pair.first, pair.second);
						if(_stored_bonded_interactions.count(pair) == 0) _stored_bonded_interactions[pair] = e;
					}
				}

				CONFIG_INFO->notify(CONFIG_INFO->box->UPDATE_EVENT);
			}
			else {
				// volume move rejected
				_box->init(old_box_sides.x, old_box_sides.y, old_box_sides.z);
				for(int k = 0; k < N(); k++) {
					BaseParticle *p = _particles[k];
					//p->pos /= _box_side / old_box_side;
					p->pos = _particles_old[k]->pos;
					p->set_ext_potential(current_step(), _box.get());
				}
				_interaction->set_is_infinite(false);
				//for (int i = 0; i < _N; i ++) _lists->single_update(_particles[i]);
				if(!_lists->is_updated()) {
					_timer_lists->resume();
					//printf ("updating lists after  rejection\n");
					_lists->global_update();
					_N_updates++;
					_timer_lists->pause();
				}
				if((current_step() < _MC_equilibration_steps && _adjust_moves) || _target_box > 0.f) _delta[MC_MOVE_VOLUME] /= 1.01;
			}
			_timer_box->pause();
		}
		else {
			_timer_move->resume();
			// do normal move
			int pi = (int) (drand48() * N());
			BaseParticle *p = _particles[pi];

			int move = (drand48() < (number) 0.5f) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;
			if(!p->is_rigid_body()) move = MC_MOVE_TRANSLATION;

			_tries[move]++;
			number delta_E = -_particle_energy(p, true);
			//number delta_E = -_particle_energy(p, false);
			p->set_ext_potential(current_step(), _box.get());
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
			_lists->single_update(p);

			//get_time(&_timer, 2);
			if(!_lists->is_updated()) {
				//printf ("updating lists because of translations\n");
				_timer_lists->resume();
				_lists->global_update();
				_N_updates++;
				_timer_lists->pause();
			}
			//get_time(&_timer, 3);

			_stored_bonded_tmp.clear();
			delta_E += _particle_energy(p, false);
			p->set_ext_potential(current_step(), _box.get());
			delta_E_ext += p->ext_potential;

			// uncomment to check the energy at a given time step.
			// may be useful for debugging purposes
			//if (curr_step > 410000 && curr_step <= 420001)
			// printf("delta_E: %lf\n", (double)delta_E);

			if(!_overlap && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / _T) > drand48())) {
				_accepted[move]++;
				_U += delta_E;
				if(current_step() < _MC_equilibration_steps && _adjust_moves) {
					_delta[move] *= 1.03;
					if(move == MC_MOVE_TRANSLATION && _delta[move] > _verlet_skin * sqrt(3.) / 2.) _delta[move] = _verlet_skin * sqrt(3.) / 2.;
					if(move == MC_MOVE_ROTATION && _delta[move] > M_PI / 2.) _delta[move] = M_PI / 2.;
				}

				// slightly faster than doing a loop over the indexes of affected
				for(auto &pair : _stored_bonded_tmp) {
					_stored_bonded_interactions[pair.first] = pair.second;
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
				_lists->single_update(p);
				_interaction->set_is_infinite(false);
				if(current_step() < _MC_equilibration_steps && _adjust_moves) _delta[move] /= 1.01;

			}
			_overlap = false;
			_timer_move->pause();
		}
	}

	if(_target_box > 0.) {
		if(fabs(_target_box / _box->box_sides()[0] - 1.) < _box_tolerance) {
			SimManager::stop = true;
			OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend) Box adjusted to %g (target %g, relative error %g)", _box->box_sides()[0], _target_box, fabs(_target_box / _box->box_sides()[0] - 1.));
		}
	}

	_mytimer->pause();
}

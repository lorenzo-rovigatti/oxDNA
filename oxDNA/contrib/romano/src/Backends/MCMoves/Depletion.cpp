/**
 * @file    Depletion.cpp
 * @date    30/ott/2018
 * @author  flavio
 *
 */
#include "Depletion.h"

#include <random>  // for poissonian distribution

template<typename number>
Depletion<number>::Depletion (){
	_delta_trs = (number) -1.f;
	_delta_rot = (number) -1.f;
	_delta_swm = (number) -1.f;
	_delta_trs_max = (number) -1.f;
	_delta_rot_max = (number) -1.f;
	_delta_swm_max = (number) -1.f;
	_ntries = -1;
	_sigma_dep = 0.5f;
	_rho_dep = -1.f;
	_tryvolume = -1.f;
	_mu_gas = 1.;
	_avg_dn = 0.;
	_n_dn = 0;
	_n_dn_n0 = 0;
	_z = (number) 0.f;
	_length = 10.f;
	_generator = std::default_random_engine();
}

template<typename number>
Depletion<number>::~Depletion () {

}

double clength (double r, double d) {
	return r*r*acos(d/r) - d*sqrt(r*r - d*d);
}

double coverlap (double r1, double r2, double d) {
	if (r2 > r1) {
		double tmp = r1;
		r1 = r2;
		r2 = tmp;
	}
	double d1 = (d*d + r1*r1 - r2*r2) / (2.*d);
	double d2 = (d*d - r1*r1 + r2*r2) / (2.*d);
	return clength(r1,d1) + clength(r2,d2);
}

template<typename number>
void Depletion<number>::init () {
	BaseMove<number>::init();
	_tryvolume = (10. + 2. * _sigma_dep) * pow((0.5 + _sigma_dep), 2) * M_PI;

	_z = exp(_mu_gas / this->_T);

	_poisson = std::poisson_distribution<int> (_z*_tryvolume);

	_ntries = _z * _tryvolume; // this will be extracted every time, this is just the average to log

	// let's compute an effective temperature
	double ovolume = coverlap (0.5 + _sigma_dep, 0.5 + _sigma_dep, 1.) - 2. * coverlap(0.5+_sigma_dep, 0.5, 1.);
	ovolume *= (_length + 2.*_sigma_dep);
	double kteq = _z*ovolume;

	if (this->_restrict_to_type < 0) throw oxDNAException("Depletion move MUST be restricted to a type");
	OX_LOG(Logger::LOG_INFO, "(Depletion.cpp) Depletion move initiated with delta_trs=%g, delta_trs_max=%g, delta_rot=%g, delta_rot_max=%g, delta_swm=%g, delta_swm_max=%g", _delta_trs, _delta_trs_max, _delta_rot, _delta_rot_max, _delta_swm, _delta_swm_max);
	OX_LOG(Logger::LOG_INFO, "(Depletion.cpp)                               <ntries>=%d, sigma_dep=%g, mu_gas=%g, tryvolume=%g", _ntries, _sigma_dep, _mu_gas, _tryvolume);
	OX_LOG(Logger::LOG_INFO, "(Depletion.cpp)                               restrict_to_type=%d, and probability %g", this->_restrict_to_type, this->prob);
	OX_LOG(Logger::LOG_INFO, "(Depletion.cpp)                               max overlap volume=%g, effective T=%g", ovolume, 1./kteq);
}

template<typename number>
void Depletion<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	double tmpf[3];
	std::string tmpstr;
	getInputString (&inp, "deltas", tmpstr, 1);
	int tmpi = sscanf(tmpstr.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	_delta_trs = (number) tmpf[0];
	_delta_rot = (number) tmpf[1];
	_delta_swm = (number) tmpf[2];
	if(tmpi != 3) throw oxDNAException ("(Depletion.cpp) Could not parse deltas (found deltas=%s, provide deltas=<float>,<float>,<float>)", tmpstr.c_str());
	getInputString (&inp, "deltas_max", tmpstr, 1);
	tmpi = sscanf(tmpstr.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException ("(Depletion.cpp) Could not parse deltas (found deltas_max=%s, provide deltas_max=<float>,<float>,<float>)", tmpstr.c_str());
	_delta_trs_max = (number) tmpf[0];
	_delta_rot_max = (number) tmpf[1];
	_delta_swm_max = (number) tmpf[2];
	getInputNumber (&inp, "sigma_dep", &_sigma_dep, 1);
	//getInputNumber (&inp, "rho_dep", &_rho_dep, 1);
	getInputNumber (&inp, "mu_gas", &_mu_gas, 1);
	getInputInt (&inp, "ntries", &_ntries, 1);
	getInputNumber (&inp, "length", &_length, 0);
}

int part_log(int n, int k) {
	if (n < k) throw oxDNAException ("never call me like that, you useless programmer");
	if (n == k) return 1;
	if (n == k + 1) return n;
	else return n * part_log(n - 1, k);
}

template<typename number>
void Depletion<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * p = this->_Info->particles[pi];
	if (this->_restrict_to_type >= 0) {
		while(p->type != this->_restrict_to_type) {
			pi = (int) (drand48() * (*this->_Info->N));
			p = this->_Info->particles[pi];
		}
	}

	// compute the energy before the move
	number delta_E;
	if (this->_compute_energy_before) delta_E = -this->particle_energy(p);
	else delta_E = (number) 0.f;
	p->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential;

	std::vector<BaseParticle<number> *> neighs_old = this->_Info->lists->get_complete_neigh_list(p);

	// move_particle
	int move_type;
	int choice = lrand48() % 21;
	if ( 0 <= choice && choice < 10) move_type = DEPLETION_TRANSLATION;
	if (10 <= choice && choice < 20) move_type = DEPLETION_ROTATION;
	if (20 <= choice && choice < 21) move_type = DEPLETION_SWIM;

	LR_vector<number> pos_old = p->pos;
	LR_matrix<number> orientation_old = p->orientation;
	LR_matrix<number> orientationT_old = p->orientationT;
	if (move_type == DEPLETION_TRANSLATION) {
		number dx = (number) 2.f * _delta_trs * (drand48() - 0.5);
		number dy = (number) 2.f * _delta_trs * (drand48() - 0.5);
		p->pos += p->orientation.v1 * dx + p->orientation.v2 * dy;
		//p->pos.z += (number) 2.f * _delta_trs * (drand48() - 0.5); // better to have separate on-axis and off-axis moves (SWIM)
	}
	if (move_type == DEPLETION_ROTATION) {
		LR_matrix<number> R = Utils::get_random_rotation_matrix_from_angle<number> (_delta_rot * drand48());
		p->orientation = R * p->orientation;
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();
	}
	if (move_type == DEPLETION_SWIM) {
		p->pos += ((number) 2.f * _delta_swm * (number)(drand48() - 0.5)) * p->orientation.v3;
	}

	this->_Info->lists->single_update(p);
	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	// see if there is an overlap
	std::vector<BaseParticle<number> *> neighs_new = this->_Info->lists->get_complete_neigh_list(p);
	typename std::vector<BaseParticle<number> *>::iterator it;

	// energy after the move
	delta_E += this->particle_energy(p);
	p->set_ext_potential(curr_step, this->_Info->box);
	delta_E_ext += p->ext_potential;

	BaseParticle<number> * p_new = p;
	BaseParticle<number> * p_old = new BaseParticle<number> ();
	p_old->type = p->type;
	p_old->pos = pos_old;
	p_old->orientation = orientation_old;
	p_old->orientationT = orientationT_old;
	p_old->index = p->index;

	// helper particle
	BaseParticle<number> * q = new BaseParticle<number> ();
	q->type = this->_restrict_to_type + 1;
	q->index = (*this->_Info->N);

	//printf ("%d N\n", *this->_Info->N);
	//throw oxDNAException ("cacca\n");

	//printf ("%d\n", _ntries);

	// if there is no overlap between large particles, we handle the depletant
	bool depletion_accept = true;
	if (this->_Info->interaction->get_is_infinite() == false) {
		// for each depletant we test if it overlaps
		//int n_dep = _rho_dep * _tryvolume; // put poissonian here

		_ntries = _poisson(_generator);

		// free volume before the move;
		int ntry = 0, nfv = 0;

//#####pragma
		while (ntry < _ntries) {
			number dx = 2.*drand48() - 1.; // between -1 and 1
			number dy = 2.*drand48() - 1.; // between -1 and 1
			number dz = drand48() - 0.5;   // between -0.5 and 0.5;
			while (dx*dx + dy*dy >= 1.) {
				dx = 2.*drand48() - 1.;
				dy = 2.*drand48() - 1.;
			}
			dx = dx * (0.5 + _sigma_dep);
			dy = dy * (0.5 + _sigma_dep);
			dz = dz * (_length + 2. * _sigma_dep);

			q->pos = p_new->pos + (p_new->orientation.v1 * dx) +
			                      (p_new->orientation.v2 * dy) +
								  (p_new->orientation.v3 * dz);

			// compare with core
			if (_sigma_dep < 0.5f) {
				if (fabs(dz) > _length / 2. - _sigma_dep) {
					if (dx * dx + dy * dy < (0.5 - _sigma_dep) * (0.5 - _sigma_dep)) {
						this->_Info->interaction->set_is_infinite(true); // to skip cycle
					}
				}
			}

			// we check if we overlap with any of the other particles, p_old included
			number junk = this->_Info->interaction->pair_interaction(p_old, q);
			for (it = neighs_new.begin(); it != neighs_new.end() && this->_Info->interaction->get_is_infinite() == false; ++it) {
				junk += this->_Info->interaction->pair_interaction(*it, q);
				if ((*it)->type != this->_restrict_to_type) throw oxDNAException("disaster %d", __LINE__);
			}

			if (this->_Info->interaction->get_is_infinite() == false) nfv ++;
			this->_Info->interaction->set_is_infinite(false);

			ntry ++;
		}

		int n_new = nfv;

		int dn = nfv;

		int avg_dn = -nfv;
		number fv_old = _tryvolume * nfv / (number) _ntries;

		ntry = 0;
		nfv = 0;
		while (ntry < _ntries) {
			number dx = 2.*drand48() - 1.; // between -1 and 1
			number dy = 2.*drand48() - 1.; // between -1 and 1
			number dz = drand48() - 0.5;   // between -0.5 and 0.5;
			while (dx*dx + dy*dy >= 1.) {
				dx = 2.*drand48() - 1.;
				dy = 2.*drand48() - 1.;
			}
			dx = dx * (0.5 + _sigma_dep);
			dy = dy * (0.5 + _sigma_dep);
			dz = dz * (_length + 2. * _sigma_dep);

			//q->pos.x = p_old->pos.x + p_old->orientation.v1.x * dx + p_old->orientation.v2.x * dy + p_old->orientation.v3.x * dz;
			//q->pos.y = p_old->pos.y + p_old->orientation.v1.y * dx + p_old->orientation.v2.y * dy + p_old->orientation.v3.y * dz;
			//q->pos.z = p_old->pos.z + p_old->orientation.v1.z * dx + p_old->orientation.v2.z * dy + p_old->orientation.v3.z * dz;

			q->pos = p_old->pos + (p_old->orientation.v1 * dx) +
								  (p_old->orientation.v2 * dy) +
								  (p_old->orientation.v3 * dz);

			// compare with core
  			if (_sigma_dep < 0.5f) {
  				if (fabs(dz) > _length / 2. - _sigma_dep) {
  					if (dx * dx + dy * dy < (0.5 - _sigma_dep) * (0.5 - _sigma_dep)) {
  						this->_Info->interaction->set_is_infinite(true); // to skip cycle
  					}
  				}
  			}

			// we check if we overlap with any of the other particles, p_old included
			number junk = this->_Info->interaction->pair_interaction(p_new, q);
			for (it = neighs_old.begin(); it != neighs_old.end() && this->_Info->interaction->get_is_infinite() == false; ++it) {
				junk += this->_Info->interaction->pair_interaction(*it, q);
				if ((*it)->type != this->_restrict_to_type) throw oxDNAException("disaster %d", __LINE__);
			}

			if (this->_Info->interaction->get_is_infinite() == false) nfv ++;
			this->_Info->interaction->set_is_infinite(false);

			ntry ++;
		}

		int n_old = nfv;

		dn = n_new - n_old;

		avg_dn += nfv;
		number fv_new = _tryvolume * nfv / (number) _ntries;

		delta_E += - _mu_gas * (fv_new  - fv_old);
		delta_E = 0.;

		//number f = pow(_z * _tryvolume, dn) * exp(lgamma(nfv + 1) - lgamma(nfv + dn + 1));
		//number f = pow(_z * _tryvolume, dn) / (number) part_log(nfv + dn, nfv);
		number f;
		if (dn >= 0) f = pow(_z * _tryvolume, dn) / part_log(n_old + dn, n_old);
		else f = pow(_z * _tryvolume, dn) * part_log(n_old, n_old + dn);

		f = 1./f;

		if (f >= (number) 1.f || f > drand48()) depletion_accept = true;
		else depletion_accept = false;

		//if (dn >= 0) printf ("dn, f: %5d %8.6lf (new=%d, old=%d, ntries=%d, part_log=%g, pow=%g)\n", dn, (double) f, n_new, n_old, _ntries, 1.0 / part_log(n_old + dn, n_old), _z * _tryvolume);
		//else printf ("dn, f: %5d %8.6lf (new=%d, old=%d, ntries=%d, part_log=%g, pow=%g)\n", dn, (double) f, n_new, n_old, _ntries, (double) part_log(n_old, n_old + dn), _z * _tryvolume);

		/*
		if (avg_dn != 0) _n_dn_n0 ++;
		_avg_dn += fabs(avg_dn);
		_n_dn ++;

		if (_n_dn == 5000) {
			number f = _avg_dn / (double) _n_dn;
			//if (f < 2) _ntries ++;
			//if (f > 4) _ntries --;
			printf ("%5.3f (%5.3f) %5d %5d (%5d)\n", _avg_dn / (double) _n_dn, _avg_dn / (double)_n_dn_n0, _n_dn, _n_dn_n0, _ntries);
			_avg_dn = 0;
			_n_dn = 0;
			_n_dn_n0 = 0;
		}*/
		//if (fabs(delta_E) > 0) printf ("%g %g\n", delta_E, delta_E / this->_T);
		//printf ("%g %g %g %g %d %d\n", delta_E, fv_new - fv_old, fv_new, fv_old, nfv, nfvold);
	}

	delete q;
	delete p_old;

	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false &&
		depletion_accept == true &&
		((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48()) ) {
		// move accepted
        this->_accepted ++;

		//if (delta_E + delta_E_ext > 0) printf ("accepted %g %g\n", delta_E, delta_E_ext);

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			if (move_type == DEPLETION_TRANSLATION) {
				_delta_trs *= this->_acc_fact;
				if (_delta_trs >= _delta_trs_max) _delta_trs = _delta_trs_max;
			}
			if (move_type == DEPLETION_ROTATION) {
				_delta_rot *= this->_acc_fact;
				if (_delta_rot >= _delta_rot_max) _delta_rot = _delta_rot_max;
			}
			if (move_type == DEPLETION_SWIM) {
				_delta_swm *= this->_acc_fact;
				if (_delta_swm >= _delta_swm_max) _delta_swm = _delta_swm_max;
			}
		}
	}
	else {
		if (move_type == DEPLETION_TRANSLATION || move_type == DEPLETION_SWIM) {
			p->pos = pos_old;
		}
		else {
			p->orientation = orientation_old;
			p->orientationT = orientationT_old;
			p->set_positions();
		}
		this->_Info->interaction->set_is_infinite(false);

		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			if (move_type == DEPLETION_TRANSLATION) _delta_trs /= this->_rej_fact;
			if (move_type == DEPLETION_ROTATION) _delta_rot /= this->_rej_fact;
			if (move_type == DEPLETION_SWIM) _delta_rot /= this->_rej_fact;
		}
	}

	return;
}

template<typename number>
void Depletion<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta_trs = %g, delta_rot = %g, delta_swm = %g", _delta_trs, _delta_rot, _delta_swm);
}

template class Depletion<float>;
template class Depletion<double>;

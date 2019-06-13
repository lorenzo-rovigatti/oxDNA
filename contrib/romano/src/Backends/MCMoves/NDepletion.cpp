/**
 * @file    NDepletion.cpp
 * @date    30/ott/2018
 * @author  flavio
 *
 */
#include "NDepletion.h"

#include <random>  // for poissonian distribution

template<typename number>
NDepletion<number>::NDepletion (){
	_ntries = -1;
	_sigma_dep = 0.5f;
	_tryvolume = -1.f;
	_mu_dep = 0.;
	_z = (number) 0.f;
	_generator = std::default_random_engine();
}

template<typename number>
NDepletion<number>::~NDepletion () {

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
void NDepletion<number>::init () {
	BaseMove<number>::init();
	_tryvolume = (10. + 2. * _sigma_dep) * pow((0.5 + _sigma_dep), 2) * M_PI;

	_z = exp(_mu_dep / this->_T);

	_poisson = std::poisson_distribution<int> (_z * _tryvolume);

	_ntries = _z * _tryvolume; // this will be extracted every time, this is just the average to log

	// let's compute an effective temperature
	double ovolume = coverlap (0.5 + _sigma_dep, 0.5 + _sigma_dep, 1.) - 2. * coverlap(0.5+_sigma_dep, 0.5, 1.);
	ovolume *= (10. + 2.*_sigma_dep);
	double kteq = _z*ovolume;

	if (this->_restrict_to_type < 0) throw oxDNAException("NDepletion move MUST be restricted to a type");
	//OX_LOG(Logger::LOG_INFO, "(NDepletion.cpp) NDepletion move initiated with delta_trs=%g, delta_trs_max=%g, delta_rot=%g, delta_rot_max=%g, delta_swm=%g, delta_swm_max=%g", _delta_trs, _delta_trs_max, _delta_rot, _delta_rot_max, _delta_swm, _delta_swm_max);
	OX_LOG(Logger::LOG_INFO, "(NDepletion.cpp) NDepletion                    <ntries>=%d, sigma_dep=%g, mu_dep=%g, tryvolume=%g", _ntries, _sigma_dep, _mu_dep, _tryvolume);
	OX_LOG(Logger::LOG_INFO, "(NDepletion.cpp)                               restrict_to_type=%d, and probability %g", this->_restrict_to_type, this->prob);
	OX_LOG(Logger::LOG_INFO, "(NDepletion.cpp)                               max overlap volume=%g, effective T=%g", ovolume, 1./kteq);
}

template<typename number>
void NDepletion<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "sigma_dep", &_sigma_dep, 1);
	getInputNumber (&inp, "mu_dep", &_mu_dep, 1);
}

int part_log(int n, int k) {
	if (n < k) throw oxDNAException ("never call me like that, you useless programmer");
	if (n == k) return 1;
	if (n == k + 1) return n;
	else return n * part_log(n - 1, k);
}

template<typename number>
void NDepletion<number>::apply (llint curr_step) {

	number length = 10.;

	// we increase the attempted count
	this->_attempted += 1;

	// we select the target particle
	int pti = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * pt = this->_Info->particles[pti];
	if (this->_restrict_to_type >= 0) {
		while(pt->type != this->_restrict_to_type) {
			pti = (int) (drand48() * (*this->_Info->N));
			pt = this->_Info->particles[pti];
		}
	}
	std::vector<BaseParticle<number> *> target_neighs = this->_Info->lists->get_complete_neigh_list(pt);

	// remove particles that are not in the bonding volume
	//typename std::vector<BaseParticle<number> *>::iterator it = target_neighs.begin();
	//while (it != target_neighs.end()) {
	//	if (!are_close(pt, *it, true)) it = target_neighs.erase(it);
	//	else ++it;
	//}

	// we select the particle to move
	int pi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * p = this->_Info->particles[pi];

	if (target_neighs.size() == 0) {
		// we may return here, since there is noone to send away
		return;
	}

	if (target_neighs.size() == 1) {
		p = target_neighs[0];
		pi = p->index;
	}
	else {
		pi = (int) (drand48() * (target_neighs.size()));
		p = target_neighs[pi];
		while (p == pt) {
			if (this->_restrict_to_type >= 0) {
				while(p->type != this->_restrict_to_type) {
					pi = (int) (drand48() * (target_neighs.size()));
					p = target_neighs[pi];
				}
			}
			else {
				pi = (int) (drand48() * (target_neighs.size()));
				p = target_neighs[pi];
			}
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
	int choice = lrand48() % 5; // 60, 120, 180, 240, 300

	LR_vector<number> pos_old = p->pos;
	LR_matrix<number> orientation_old = p->orientation;
	LR_matrix<number> orientationT_old = p->orientationT;

	LR_matrix<number> R;
	{
		LR_vector<number> axis = pt->orientation.v3;
		axis.normalize();
		number t = (1 + choice) * M_PI / 3.;

		number sintheta = sin(t);
		number costheta = cos(t);
		number olcos = 1. - costheta;

		number xyo = axis.x * axis.y * olcos;
		number xzo = axis.x * axis.z * olcos;
		number yzo = axis.y * axis.z * olcos;
		number xsin = axis.x * sintheta;
		number ysin = axis.y * sintheta;
		number zsin = axis.z * sintheta;

		R = LR_matrix<number> (axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);
	}

	p->pos = pt->pos + R * this->_Info->box->min_image(pt, p);
	p->orientation = R * p->orientation;
	p->orientationT = p->orientation.get_transpose();
	p->set_positions();

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

	//since lists can change, and after moving p could in principle not be
	// a neighbour of pt anymore, we need this check to ensure detailed balance
	if (target_neighs.size() != this->_Info->lists->get_complete_neigh_list(pt).size()) {
		//OX_LOG(Logger::LOG_DEBUG, "Not a neighbour anymore...");
		this->_Info->interaction->set_is_infinite(true);
	}

	// if there is no overlap between large particles, we handle the depletant
	bool depletion_accept = true;
	if (this->_Info->interaction->get_is_infinite() == false) {
		// for each depletant we test if it overlaps
		//int n_dep = _rho_dep * _tryvolume; // put poissonian here

		_ntries = _poisson(_generator);

		// free volume before the move;
		int ntry = 0, nfv = 0;
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
			dz = dz * (length + 2. * _sigma_dep);

			q->pos = p_new->pos + (p_new->orientation.v1 * dx) +
			                      (p_new->orientation.v2 * dy) +
								  (p_new->orientation.v3 * dz);

			// compare with core
			if (_sigma_dep < 0.5f) {
				if (fabs(dz) > length / 2. - _sigma_dep) {
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
			dz = dz * (length + 2. * _sigma_dep);

			q->pos = p_old->pos + (p_old->orientation.v1 * dx) +
								  (p_old->orientation.v2 * dy) +
								  (p_old->orientation.v3 * dz);

			// compare with core
  			if (_sigma_dep < 0.5f) {
  				if (fabs(dz) > length / 2. - _sigma_dep) {
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

		//number f = pow(_z * _tryvolume, dn) * exp(lgamma(nfv + 1) - lgamma(nfv + dn + 1));
		//number f = pow(_z * _tryvolume, dn) / (number) part_log(nfv + dn, nfv);
		number f;
		if (dn >= 0) f = pow(_z * _tryvolume, dn) / part_log(n_old + dn, n_old);
		else f = pow(_z * _tryvolume, dn) * part_log(n_old, n_old + dn);

		f = 1./f;

		if (f >= (number) 1.f || f > drand48()) depletion_accept = true;
		else depletion_accept = false;
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
		/*
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
		}*/
	}
	else {
		p->pos = pos_old;
		p->orientation = orientation_old;
		p->orientationT = orientationT_old;
		p->set_positions();
		this->_Info->interaction->set_is_infinite(false);

		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}

		/*
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			if (move_type == DEPLETION_TRANSLATION) _delta_trs /= this->_rej_fact;
			if (move_type == DEPLETION_ROTATION) _delta_rot /= this->_rej_fact;
			if (move_type == DEPLETION_SWIM) _delta_rot /= this->_rej_fact;
		}*/
	}

	return;
}

template<typename number>
void NDepletion<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	//OX_LOG(Logger::LOG_INFO, "\tdelta_trs = %g, delta_rot = %g, delta_swm = %g", _delta_trs, _delta_rot, _delta_swm);
}

template class NDepletion<float>;
template class NDepletion<double>;

/**
 * @file    DepletionVolume.cpp
 * @date    30/ott/2018
 * @author  flavio
 *
 */
#include "DepletionVolume.h"

#include <random>  // for poissonian distribution

template<typename number>
DepletionVolume<number>::DepletionVolume (){
	_delta = (number) -1.f;
	_delta_max = (number) -1.f;
	_ntries = -1;
	_sigma_dep = 0.5f;
	_mu_gas = 1.;
	_z = (number) 0.f;
	_length = 10.f;
	_pressure = -1.f;
	_change_volume = false;
	_affected_side = -1; 
	_generator = std::default_random_engine();
}

template<typename number>
DepletionVolume<number>::~DepletionVolume () {

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
void DepletionVolume<number>::init () {
	BaseMove<number>::init();

	if (_pressure < 0.f) // defaults to -1
		_change_volume = false;
	else
		_change_volume = true;

	_pos_old.resize(*this->_Info->N);
	if (_delta_max < (number) 0.) _delta_max = 1.e6;
	_z = exp(_mu_gas / this->_T);

	// let's compute an effective temperature
	double ovolume = coverlap (0.5 + _sigma_dep, 0.5 + _sigma_dep, 1.) - 2. * coverlap(0.5+_sigma_dep, 0.5, 1.);
	ovolume *= (_length + 2.*_sigma_dep);
	double kteq = _z*ovolume;


	// take a look at the affected sides of the box.
	if (_affected_side < -1 || _affected_side > 2)
		throw oxDNAException ("(DepletionVolume.cpp) ERROR: affected_side must be -1(all, default behaviour), 0 (x), 1(y) or 2(z); got wrong value %d", _affected_side);

	OX_LOG(Logger::LOG_INFO, "(DepletionVolume.cpp) DepletionVolume move initiated with pressure=%g, affected_side=%d, delta=%g, delta_max=%g", _pressure, _affected_side, _delta, _delta_max);
	OX_LOG(Logger::LOG_INFO, "(DepletionVolume.cpp)                               sigma_dep=%g, mu_gas=%g", _sigma_dep, _mu_gas);
	OX_LOG(Logger::LOG_INFO, "(DepletionVolume.cpp)                               max overlap volume=%g, effective T=%g", ovolume, 1./kteq);
}

template<typename number>
void DepletionVolume<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "delta_max", &_delta_max, 0);
	getInputNumber (&inp, "sigma_dep", &_sigma_dep, 1);
	getInputNumber (&inp, "mu_gas", &_mu_gas, 1);
	getInputNumber (&inp, "length", &_length, 0);
	getInputNumber (&inp, "pressure", &_pressure, 0);
	getInputInt (&inp, "affected_side", &_affected_side, 0);

}

int part_log(int n, int k) {
	if (n < k) throw oxDNAException ("(DepletionVolume.cpp) never call me like that, you useless programmer");
	if (n == k) return 1;
	if (n == k + 1) return n;
	else return n * part_log(n - 1, k);
}

double part_log_n(int n, int k) {
	if (n < k) return (double) part_log(k, n);
	else return 1. / (double) part_log(n, k);
}

template<typename number>
void DepletionVolume<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	BaseParticle<number> ** particles = this->_Info->particles;
	int N = *(this->_Info->N);

	LR_vector<number> box_sides = this->_Info->box->box_sides();
	LR_vector<number> old_box_sides = box_sides;

	number oldE;
	if (this->_compute_energy_before) oldE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	else oldE = (number) 0.f;
	if (this->_Info->interaction->get_is_infinite()) printf ("WHAAT\n");

	number oldV = this->_Info->box->V();

	if (_change_volume == false) {
		number dL = _delta * (drand48() - 0.5);
		box_sides.z += dL;

		number f = sqrt(old_box_sides.z / box_sides.z);
		box_sides.y = f * box_sides.y;
		box_sides.x = f * box_sides.x;
	}
	else {
		// random walk in log(V)
		number dv = 2. * _delta * (drand48() - 0.5);
		if (_affected_side == -1) {
			// all sides affected in the same way
			number f = exp(dv/3.);
			box_sides.x *= f;
			box_sides.y *= f;
			box_sides.z *= f;
		}
		else {
			number f = exp(dv);
			box_sides[_affected_side] *= f;
		}
	}

	RodCells<number> * cells = static_cast<RodCells<number> * > (this->_Info->lists);
	int old_ncells = cells->get_N_cells();
	std::vector<std::vector<BaseParticle<number> *> >old_occupation;
	for (int k = 0; k < old_ncells; k ++) {
		old_occupation.push_back(cells->whos_there(k));
	}

	this->_Info->box->init(box_sides[0], box_sides[1], box_sides[2]);
	number dExt = (number) 0.f;
	for(int k = 0; k < N; k ++) {
		BaseParticle<number> *p = particles[k];
		dExt -= p->ext_potential;
		_pos_old[k] = p->pos;
		p->pos.x *= box_sides.x/old_box_sides.x;
		p->pos.y *= box_sides.y/old_box_sides.y;
		p->pos.z *= box_sides.z/old_box_sides.z;
		p->set_ext_potential(curr_step, this->_Info->box);
		dExt += p->ext_potential;
	}

	// this bit has to come after the update of particles' positions
	if(!this->_Info->lists->is_updated()) this->_Info->lists->global_update();

	int new_ncells = cells->get_N_cells();
	std::vector<std::vector<BaseParticle<number> * > > new_occupation;
	if (new_ncells != old_ncells) {
		for (int k = 0; k < new_ncells; k ++) {
			new_occupation.push_back(cells->whos_there(k));
		}
	}
	else {
		new_occupation.resize(new_ncells);
		for (int k = 0; k < new_ncells; k ++) {
			new_occupation[k] = old_occupation[k];
		}
	}

	number newE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	number dE = newE - oldE + dExt;
	number V = this->_Info->box->V();
	number dV = V - oldV; 

	if (_change_volume == false && fabs(dV/V) > 1.e-3) throw oxDNAException ("dV/V too large: %g", dV/V);

	int n = 0;
	if (this->_Info->interaction->get_is_infinite() == false && V < oldV) {
		// fake particle
		BaseParticle<number> * p = new BaseParticle<number> ();
		p->orientation.v1 = LR_vector<number> (0., 0., 1.);
		p->orientation.v2 = LR_vector<number> (0., 1., 0.);
		p->orientation.v3 = LR_vector<number> (1., 0., 0.);
		p->orientationT = p->orientation.get_transpose();
		p->index = N;
		p->type = 1;

		// old conf

		// we could do this before the move, but we reject most moves because of cylinder overlaps,
		// so it is best to do it only if no overlaps occur
		_poisson = std::poisson_distribution<int> (_z * oldV);
		int n_old = 0;
		int ndep = _poisson(_generator);
		std::vector<LR_vector<number> > dep_pos;
		this->_Info->box->init(old_box_sides[0], old_box_sides[1], old_box_sides[2]);
		if(!this->_Info->lists->is_updated()) this->_Info->lists->global_update();
		for (int d = 0; d < ndep; d++) {
			p->pos.x = drand48() * old_box_sides.x;
			p->pos.y = drand48() * old_box_sides.y;
			p->pos.z = drand48() * old_box_sides.z;
			bool overlap = false;
			int k = cells->get_cell_index(p->pos);
			typename std::vector<BaseParticle<number> *>::iterator it;
			for (it = old_occupation[k].begin(); it != old_occupation[k].end() && overlap == false; ++it) {
				BaseParticle<number> * q = *it;
				LR_vector<number> dr = this->_Info->box->min_image(_pos_old[q->index], p->pos);
				number rdot = dr * q->orientation.v3;
				if (fabs(rdot) <= 0.5 * _length) {
					LR_vector<number> r_on_axis = rdot * q->orientation.v3;
					LR_vector<number> r_off_axis = dr - r_on_axis;
					if (r_off_axis.norm() <= (0.5+_sigma_dep)*(0.5+_sigma_dep)) {
						overlap = true;
					}
				}
			}
			if (overlap == false) {
				dep_pos.push_back(p->pos);
				n_old ++;
			}
		}

		// new conf
		this->_Info->box->init(box_sides[0], box_sides[1], box_sides[2]);
		if(!this->_Info->lists->is_updated()) this->_Info->lists->global_update();
		number fact = box_sides.x / old_box_sides.x;
		for (int d = 0; d < n_old; d++) {
			p->pos =  fact * dep_pos[d];
			int k = cells->get_cell_index(p->pos);
			bool overlap = false;
			typename std::vector<BaseParticle<number> *>::iterator it;
			for (it = new_occupation[k].begin(); it != new_occupation[k].end() && overlap == false; ++it) {
				BaseParticle<number> * q = *it;
				LR_vector<number> dr = this->_Info->box->min_image(q, p);  // q is the rod
				number rdot = dr * q->orientation.v3;
				if (fabs(rdot) <= 0.5 * _length) {
					LR_vector<number> r_on_axis = rdot * q->orientation.v3;
					LR_vector<number> r_off_axis = dr - r_on_axis;
					if (r_off_axis.norm() <= (0.5+_sigma_dep)*(0.5+_sigma_dep)) {
						overlap = true;
					}
				}
			}
			if (overlap == true) {
				this->_Info->interaction->set_is_infinite(true);
				break;
			}
		}

		n = n_old;
		delete p;

		//dE += -dN * log(_z) - n_new * log(V) + n_old * log(oldV) - lgamma(n_old + 1) + lgamma(n_new + 1);
	}

	// accept or reject?
	// for the ideal gas, beta * P = rho = z
	//exp(-(dE + _P*dV - (N + 1)*this->_T*log(V/oldV))/this->_T
	if (this->_Info->interaction->get_is_infinite() == false &&
		exp(-(dE + dExt + _pressure * dV) / this->_T + (N + n + 1)*log(V/oldV)) > drand48() ) {
		// move accepted
        this->_accepted ++;

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta *= this->_acc_fact;
			if (_delta >= _delta_max) _delta = _delta_max;
		}
	}
	else {
		for (int k = 0; k < N; k ++) {
			BaseParticle<number> * p = particles[k];
			p->pos = _pos_old[k];
			p->set_ext_potential(curr_step, this->_Info->box);
		}
		this->_Info->interaction->set_is_infinite(false);
		this->_Info->box->init(old_box_sides[0], old_box_sides[1], old_box_sides[2]);
		if(!this->_Info->lists->is_updated()) this->_Info->lists->global_update();

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta /= this->_rej_fact;
		}
	}

	return;
}

template<typename number>
void DepletionVolume<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}

template class DepletionVolume<float>;
template class DepletionVolume<double>;

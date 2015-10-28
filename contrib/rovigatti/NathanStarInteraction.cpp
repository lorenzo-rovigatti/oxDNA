/*
 * NathanStarInteraction.cpp
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#include "NathanStarInteraction.h"

template<typename number>
NathanStarInteraction<number>::NathanStarInteraction() : _sqrt_pi(sqrt(M_PI)), _patch_r(0.5), _sqr_patch_r(SQR(_patch_r)), _star_f(-1), _is_marzi(false), _make_crystal(false), _interp_size(500) {
	this->_int_map[PATCHY_PATCHY] = &NathanStarInteraction<number>::_patchy_interaction;

	_rep_E_cut = 0.;
	_rep_power = 200;

	_patch_alpha = 0.12;
	_patch_power = 30;
	_star_factor = 1.;

	_N_stars = _N_patchy = 0;

	_spl_patchy = NULL;
	_acc_patchy = NULL;
	_spl_patchy_star = NULL;
	_acc_patchy_star = NULL;
}

template<typename number>
NathanStarInteraction<number>::~NathanStarInteraction() {
	if(_spl_patchy_star != NULL) {
		gsl_spline_free(_spl_patchy);
		gsl_interp_accel_free(_acc_patchy);

		gsl_spline_free(_spl_patchy_star);
		gsl_interp_accel_free(_acc_patchy_star);
	}
}

template<typename number>
int NathanStarInteraction<number>::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d\n", &_N_patchy, &_N_stars);
	return _N_patchy + _N_stars;
}

template<typename number>
void NathanStarInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	// the number of "strands" is given by the number of chains + the number of patchy particles
	// since those are not linked to anything else
	*N_strands = _N_stars + _N_patchy;
	allocate_particles(particles, N);
}

template<typename number>
void NathanStarInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < _N_patchy; i++) {
		NathanPatchyParticle<number> *new_p = new NathanPatchyParticle<number>();
		new_p->index = i;
		new_p->type = new_p->btype = PATCHY_PARTICLE;
		new_p->strand_id = i;

		particles[i] = new_p;
	}

	for(int i = _N_patchy; i < N; i++) {
		NathanPolymerParticle<number> *new_p = new NathanPolymerParticle<number>();
		new_p->index = i;
		new_p->type = new_p->btype = POLYMER;
		new_p->strand_id = i;
		particles[i] = new_p;
	}
}

template<typename number>
void NathanStarInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "T", &_T, 1);

	getInputNumber(&inp, "NATHAN_alpha", &_patch_alpha, 0);
	getInputNumber(&inp, "NATHAN_cosmax", &_patch_cosmax, 1);

	getInputInt(&inp, "NATHAN_interp_size", &_interp_size, 0);

	number size_ratio;
	getInputNumber(&inp, "NATHAN_size_ratio", &size_ratio, 1);
	_star_sigma_g = 1. / size_ratio;

	if(_star_sigma_g > 1.) _is_marzi = true;

	getInputInt(&inp, "NATHAN_f", &_star_f, 1);
	getInputNumber(&inp, "NATHAN_star_factor", &_star_factor, 0);

	getInputBool(&inp, "NATHAN_make_crystal", &_make_crystal, 0);
	if(_make_crystal) {
		getInputString(&inp, "NATHAN_crystal_type", _crystal_type, 1);
		_N_in_crystal = -1;
		getInputInt(&inp, "NATHAN_N_in_crystal", &_N_in_crystal, 0);
	}
}

template<typename number>
void NathanStarInteraction<number>::_setup_lambda_kappa() {
	if(_star_f < 2 || _star_f > 1e5) throw oxDNAException("NATHAN_f (%d) should be > 1 and < 10e5", _star_f);

	int fs[12] = {2, 5, 10, 15, 18, 30, 40, 50, 65, 80, 100, 100000};
	number lambdas[12] = {0.46, 0.35, 0.30, 0.28, 0.27, 0.24, 0.24, 0.23, 0.23, 0.22, 0.22, 0.139};
	number kappas[12] = {0.58, 0.68, 0.74, 0.76, 0.77, 0.83, 0.85, 0.86, 0.87, 0.88, 0.89, 1.};

	int idx = 0;
	while(fs[idx] < _star_f) idx++;

	if(fs[idx] == _star_f) {
		_pi_lambda = lambdas[idx];
		_kappa_rs = kappas[idx];
	}
	else {
		number delta_x = (fs[idx] - fs[idx - 1]);
		number delta_f = (_star_f - fs[idx - 1]);
		number delta_l = lambdas[idx] - lambdas[idx - 1];
		number delta_k = kappas[idx] - kappas[idx - 1];

		_pi_lambda = lambdas[idx-1] + delta_l/delta_x*delta_f;
		_kappa_rs = kappas[idx-1] + delta_k/delta_x*delta_f;
	}
}

template<typename number>
number NathanStarInteraction<number>::_psi_1(number x, number smax) {
	return _xi * _exp_sqr_kappa_rs * (exp(-SQR(_kappa*x))/x - exp(-SQR(_kappa*smax))/smax) / _star_rs;
}

template<typename number>
number NathanStarInteraction<number>::_psi_2(number x, number smax) {
	return _xi * _exp_sqr_kappa_rs * (_sqrt_pi*(erf(_kappa*smax) - erf(_kappa*x))/_kappa + x*exp(-SQR(_kappa*x)) - smax*exp(-SQR(_kappa*smax))) / _star_rs;
}

template<typename number>
number NathanStarInteraction<number>::_patchy_star_derivative(number r) {
	number sqr_r = r*r;
	number z = r - _patch_r;
	number smax = sqrt(z*(z + 1.));

	number derivative = -_T * _pi_lambda * _star_f3_2 * _patch_r / sqr_r;

	if(z < _star_rs) derivative *= (sqr_r - _sqr_patch_r)*(0.5/SQR(z) - 0.5/_sqr_star_rs + _psi_1(_star_rs, smax)) - log(z/_star_rs) + _psi_2(_star_rs, smax);
	else derivative *= (sqr_r - _sqr_patch_r)*_psi_1(z, smax) + _psi_2(z, smax);

	return derivative;
}

template<typename number>
number NathanStarInteraction<number>::_patchy_star_marzi_derivative(number r, gsl_spline *spl, gsl_interp_accel *acc, number bin) {
	number sqr_r = r*r;
	number z = r - _patch_r;
	number factor = -M_PI * _patch_r / sqr_r;

	if(spl == NULL) {
		spl = gsl_spline_alloc(gsl_interp_cspline, _interp_size);
		acc = gsl_interp_accel_alloc();

		double *ss = new double[_interp_size]();
		double *integrand = new double[_interp_size]();
		for(int i = 0; i < _interp_size; i++) {
			double s = (i + 0.5)*bin;
			double t = (z*(r + _patch_r) - SQR(s)) / s;
			ss[i] = s;
			integrand[i] = (sqr_r -_sqr_patch_r + SQR(s)) * (_pressure(s) - _pressure(s + t));
		}
		gsl_spline_init(spl, ss, integrand, _interp_size);

		delete[] ss;
		delete[] integrand;
	}

	number smax = sqrt(z*(z + 1.));
	number smax_used = (_interp_size-0.5)*bin;
	if(smax > smax_used) smax = smax_used;
	if(z > smax) return 0.f;

	return factor * gsl_spline_eval_integ(spl, z, smax, acc);
}

template<typename number>
number NathanStarInteraction<number>::_pressure(number s) {
	number common_factor = _T * _pi_lambda * _star_f3_2 / M_PI;
	if(s < _star_rs) return common_factor / CUB(s);

	number sqr_s = SQR(s);
	number sqr_kappa = SQR(_kappa);
	return common_factor * (1./sqr_s + 2*sqr_kappa) * _xi * exp(-sqr_kappa*(sqr_s - _sqr_star_rs)) / _star_rs;
}

template<typename number>
void NathanStarInteraction<number>::_setup_interp() {
	_spl_patchy = gsl_spline_alloc(gsl_interp_cspline, _interp_size);
	_acc_patchy = gsl_interp_accel_alloc();

	_spl_patchy_star = gsl_spline_alloc(gsl_interp_cspline, _interp_size);
	_acc_patchy_star = gsl_interp_accel_alloc();

	// patchy-star
	number bin = (_patchy_star_rcut - _patch_r) / _interp_size;
	number r = _patch_r + bin;
	double *rs = new double[_interp_size]();
	double *f_ys = new double[_interp_size]();
	gsl_spline *marzi_spline = NULL;
	gsl_interp_accel *marzi_accel = NULL;
	for(int i = 0; i < _interp_size; i++, r += bin) {
		rs[i] = r;
		if(_is_marzi) f_ys[i] = _patchy_star_marzi_derivative(r, marzi_spline, marzi_accel, bin);
		else f_ys[i] = _patchy_star_derivative(r);
	}
	gsl_spline_init(_spl_patchy_star, rs, f_ys, _interp_size);

	// print out the "real" derivative, its interpolation and their relative difference
	ofstream out_der("tabulated_patchy_star.dat");
	for(int i = 0; i < _interp_size-1; i++) {
		r = rs[0] + 0.001 + bin*i;
		number real_der;
		if(_is_marzi) real_der = _patchy_star_marzi_derivative(r, marzi_spline, marzi_accel, bin);
		else real_der = _patchy_star_derivative(r);
		number interp_der = gsl_spline_eval(_spl_patchy_star, r, _acc_patchy_star);
		out_der << r << " " << real_der << " " << interp_der << " " << (real_der - interp_der) / real_der << endl;
	}
	out_der.close();

	double *u_ys = new double[_interp_size]();
	double u_shift = gsl_spline_eval_integ(_spl_patchy_star, rs[0], rs[_interp_size-1], _acc_patchy_star);
	for(int i = 0; i < _interp_size; i++) {
		u_ys[i] = gsl_spline_eval_integ(_spl_patchy_star, rs[0], rs[i], _acc_patchy_star) - u_shift;
	}
	gsl_spline_init(_spl_patchy_star, rs, u_ys, _interp_size);

	delete[] rs;
	delete[] f_ys;
	delete[] u_ys;

	if(_is_marzi) {
		gsl_spline_free(marzi_spline);
		gsl_interp_accel_free(marzi_accel);
	}
}

template<typename number>
void NathanStarInteraction<number>::init() {
	_setup_lambda_kappa();

	_star_f1_2 = sqrt(_star_f);
	_star_f3_2 = CUB(_star_f1_2);

	_star_rg = 0.5*_star_sigma_g;
	_star_rs = 2.*_star_rg/3.;
	_star_sigma_s = 2.*_star_sigma_g/3.;
	_sqr_star_rs = SQR(_star_rs);
	_sqr_star_sigma_s = SQR(_star_sigma_s);
	_sqr_star_sigma_g = SQR(_star_sigma_g);

	_patchy_star_rcut = _patch_r + _star_sigma_g;
	_sqr_patchy_star_rcut = SQR(_patchy_star_rcut);

	_star_rcut = 3.*_star_sigma_s;
	_sqr_star_rcut = SQR(_star_rcut);

	_kappa = _kappa_rs / _star_rs;
	_xi = 1. / (1. + 2.*SQR(_kappa_rs));
	_exp_sqr_kappa_rs = exp(SQR(_kappa_rs));
	_zeta = _sqrt_pi * _xi * erfc(_kappa_rs) * _exp_sqr_kappa_rs / _kappa_rs;
	_erfc_kappa_rs = erfc(_kappa_rs);

	_patch_cutoff = _patch_alpha * 1.5;
	_patchy_rcut = 1. + _patch_cutoff;
	_sqr_patchy_rcut = SQR(_patchy_rcut);
	this->_rcut = 1. + _patch_cutoff;
	if(_star_rcut > this->_rcut) this->_rcut = _star_rcut;
	this->_sqr_rcut = SQR(this->_rcut);
	_patch_angular_cutoff = _patch_cosmax * 1.5;

	_rep_E_cut = pow((number) this->_rcut, -_rep_power);

	_patch_pow_sigma = pow(_patch_cosmax, _patch_power);
	_patch_pow_alpha = pow(_patch_alpha, (number) 10.);

	_setup_interp();

	OX_LOG(Logger::LOG_INFO, "pi*lambda: %lf, kappa*rs: %lf, zeta: %lf, xi: %lf", _pi_lambda, _kappa_rs, _zeta, _xi);
	OX_LOG(Logger::LOG_INFO, "patchy rcut: %lf, patchy-star rcut: %lf, star-star rcut: %lf", _patchy_rcut, _patchy_star_rcut, _star_rcut);
}

template<typename number>
void NathanStarInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
number NathanStarInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number NathanStarInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.;
}

template<typename number>
number NathanStarInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

//	p->pos = LR_vector<number>(0, 0, 0);
//	for(number mr = 0.501; mr < _patchy_star_rcut; mr += 0.001) {
//		q->pos = LR_vector<number>(mr, 0, 0);
//		computed_r = this->_box->min_image(p->pos, q->pos);
//		r = &computed_r;
//		printf("%lf %lf\n", mr, _patchy_star_interaction(p, q, r, update_forces));
//		_patchy_star_interaction(p, q, r, update_forces);
//	}
//	throw oxDNAException("fatto");

	int type = p->type + q->type;
	if(type == PATCHY_PATCHY) return _patchy_interaction(p, q, r, update_forces);
	else if(type == PATCHY_POLYMER) return _patchy_star_interaction(p, q, r, update_forces);
	else return _star_star_interaction(p, q, r, update_forces);
}

template<typename number>
void NathanStarInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N, number box_side) {
	if(!_make_crystal) {
		BaseInteraction<number, NathanStarInteraction<number> >::generate_random_configuration(particles, N, box_side);
		return;
	}
	// from here on it is implied that _make_crystal == true

	if(_N_in_crystal == -1) _N_in_crystal = _N_patchy;

	// check on the number of particles
	if(_N_in_crystal % 4 != 0) throw oxDNAException("The number of patchy particles should be a multiple of 4");
	int cxcycz = _N_in_crystal / 4;
	int cell_side = round(cbrt(cxcycz));
	if(CUB(cell_side) != cxcycz) throw oxDNAException("The number of patchy particles should be equal to 4*i^3, where i is an integer. For this specific case you may want to use %d", 4*CUB(cell_side));

	OX_LOG(Logger::LOG_INFO, "Crystal cell_side: %d", cell_side);

	number gap = 0.05;
	number sigma_c = 1.;

	LR_vector<number> tetrahedron_up[4];
	tetrahedron_up[0].x = tetrahedron_up[0].y = tetrahedron_up[0].z = 0;

	tetrahedron_up[1].x = sigma_c+gap;
	tetrahedron_up[1].y = tetrahedron_up[1].z = 0;

	tetrahedron_up[2].x = 0.5*(sigma_c+gap);
	tetrahedron_up[2].y = sqrt(3.0)/2.0*(sigma_c+gap);
	tetrahedron_up[2].z = 0;

	tetrahedron_up[3].x = 0.5*(sigma_c+gap);
	tetrahedron_up[3].y = 1./3.*sqrt(3.0)/2.0*(sigma_c+gap);
	tetrahedron_up[3].z = sqrt(2.0/3.0)*(sigma_c+gap);

	LR_matrix<number> tetrahedron_orientT[4];
	for(int i = 0; i < 4; i++) {
		LR_vector<number> axis(0, 0, 0);
		for(int j = 0; j < 4; j++) {
			if(j != i) axis += tetrahedron_up[j];
		}
		axis /= 3;
		axis -= tetrahedron_up[i];
		axis.normalize();
		tetrahedron_orientT[i] = Utils::get_random_rotation_matrix_from_angle<number>(M_PI);
		tetrahedron_orientT[i].v3 = axis;
		tetrahedron_orientT[i].orthonormalize();
	}

	int cells[3] = {cell_side, cell_side, cell_side};

	// let's build the crystal
	if(_crystal_type == "hT") {
		LR_vector<number> tetrahedron_down[4] = {tetrahedron_up[0], tetrahedron_up[1], tetrahedron_up[2], tetrahedron_up[3]};
		tetrahedron_down[2].y *= -1;
		tetrahedron_down[3].y *= -1;

		number xshift = (sigma_c + gap);
		number yshift = sqrt(3.0)*(sigma_c + gap);
		number zshift = 2.0*sqrt(2.0/3.0)*(sigma_c + gap);

		LR_matrix<number> tetrahedron_down_orientT[4];
		for(int i = 0; i < 4; i++) {
			LR_vector<number> axis(0, 0, 0);
			for(int j = 0; j < 4; j++) {
				if(j != i) axis += tetrahedron_down[j];
			}
			axis /= 3;
			axis -= tetrahedron_down[i];
			axis.normalize();
			tetrahedron_down_orientT[i] = Utils::get_random_rotation_matrix_from_angle<number>(M_PI);
			tetrahedron_down_orientT[i].v3 = axis;
			tetrahedron_down_orientT[i].orthonormalize();
		}

		int inserted = 0;
		for(int zc = 0; zc < cells[2]; zc++) {
			LR_vector<number> *tetra = (zc % 2 == 0) ? tetrahedron_up : tetrahedron_down;
			LR_matrix<number> *tetra_orientT = (zc % 2 == 0) ? tetrahedron_orientT : tetrahedron_down_orientT;
			number tot_zshift = zc*zshift;
			for(int yc = 0; yc < cells[1]; yc++) {
				number tot_yshift = yc*yshift;
				if(zc % 2 != 0) tot_yshift += yshift;
				for(int xc = 0; xc < cells[2]; xc++) {
					number tot_xshift = (yc%2 + 2*xc)*xshift;
					// loop over the tetrahedron
					for(int t = 0; t < 4; t++, inserted++) {
						BaseParticle<number> *p = particles[inserted];
						p->pos.x = tetra[t].x + tot_xshift;
						p->pos.y = tetra[t].y + tot_yshift;
						p->pos.z = tetra[t].z + tot_zshift;

						p->orientationT = tetra_orientT[t];
						p->orientation = p->orientationT.get_transpose();
						p->set_positions();
					}
				}
			}
		}
	}
	else if(_crystal_type == "cT") {
		number xshift = 2.0*(sigma_c + gap);
		number yshift = 2.0*sqrt(3.0)/2.0*(sigma_c + gap);
		number zshift = 2.0*sqrt(2.0/3.0)*(sigma_c + gap);

		int inserted = 0;
		for(int zc = 0; zc < cells[2]; zc++) {
			number tot_zshift = 3*zc*zshift;
			tot_zshift += (zc % 3)*zshift;
			tot_zshift = zc*zshift;

			for(int yc = 0; yc < cells[1]; yc++) {
				number tot_yshift = yc*yshift;
				if(zc % 3 == 1) tot_yshift -= 2.*yshift/3.;
				else if(zc % 3 == 2) tot_yshift += 2.*yshift/3.;
				for(int xc = 0; xc < cells[2]; xc++) {
					number tot_xshift = xc*xshift;
					if(yc % 2 == 1) tot_xshift += 0.5*xshift;
					// loop over the tetrahedron
					for(int t = 0; t < 4; t++, inserted++) {
						BaseParticle<number> *p = particles[inserted];
						p->pos.x = tetrahedron_up[t].x + tot_xshift;
						p->pos.y = tetrahedron_up[t].y + tot_yshift;
						p->pos.z = tetrahedron_up[t].z + tot_zshift;

						p->orientationT = tetrahedron_orientT[t];
						p->orientation = p->orientationT.get_transpose();
						p->set_positions();
						p->set_ext_potential(0, box_side);
					}
				}
			}
		}
	}
	else throw oxDNAException("The crystal type be should be either 'hT' or 'cT' and not '%s'", _crystal_type.c_str());

	// now we add the polymers
	Cells<number> c(N, this->_box);
	c.init(particles, this->_rcut);

	for(int i = _N_in_crystal; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		p->pos = LR_vector<number>(0., 0., 0.);
	}

	c.global_update();

	for(int i = _N_in_crystal; i < N; i++) {
		BaseParticle<number> *p = particles[i];

		bool inserted = false;
		do {
			p->pos = LR_vector<number> (drand48()*box_side, drand48()*box_side, drand48()*box_side);
			// random orientation
			p->orientation = Utils::get_random_rotation_matrix_from_angle<number> (M_PI);
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();

			p->set_positions();
			p->set_ext_potential(0, box_side);
			c.single_update(p);

			inserted = true;

			// we take into account the bonded neighbours
			for (unsigned int n = 0; n < p->affected.size(); n ++) {
				BaseParticle<number> * p1 = p->affected[n].first;
				BaseParticle<number> * p2 = p->affected[n].second;
				number e = 0.;
				if (p1->index <= p->index && p2->index <= p->index) {
					e = pair_interaction_bonded(p1, p2);
				}
				if(e > this->_energy_threshold) inserted = false;
			}

			// here we take into account the non-bonded interactions
			vector<BaseParticle<number> *> neighs = c.get_neigh_list(p, true);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle<number> *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < p->index && this->generate_random_configuration_overlap(p, q, box_side)) inserted = false;
			}

			// we take into account the external potential
			if (p->ext_potential > this->_energy_threshold) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && i % (N/10) == 0) OX_LOG(Logger::LOG_INFO, "Inserted %d%% of the particles (%d/%d)", i*100/N, i, N);
	}
}

template class NathanStarInteraction<float>;
template class NathanStarInteraction<double>;

extern "C" NathanStarInteraction<float> *make_interaction_float() { return new NathanStarInteraction<float>(); }
extern "C" NathanStarInteraction<double> *make_interaction_double() { return new NathanStarInteraction<double>(); }

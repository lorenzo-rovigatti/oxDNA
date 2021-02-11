/*
 * NathanStarInteraction.cpp
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#include "NathanStarInteraction.h"

NathanStarInteraction::NathanStarInteraction() :
				_sqrt_pi(sqrt(M_PI)),
				_patch_r(0.5),
				_sqr_patch_r(SQR(_patch_r)),
				_star_f(-1),
				_is_marzi(false),
				_make_crystal(false),
				_interp_size(500) {
	this->_int_map[PATCHY_PATCHY] = &NathanStarInteraction::_patchy_interaction;

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

NathanStarInteraction::~NathanStarInteraction() {
	if(_spl_patchy_star != NULL) {
		gsl_spline_free(_spl_patchy);
		gsl_interp_accel_free(_acc_patchy);

		gsl_spline_free(_spl_patchy_star);
		gsl_interp_accel_free(_acc_patchy_star);
	}
}

int NathanStarInteraction::get_N_from_topology() {
	char line[512];
	std::ifstream topology;
	topology.open(this->_topology_filename, std::ios::in);
	if(!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d\n", &_N_patchy, &_N_stars);
	return _N_patchy + _N_stars;
}

void NathanStarInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	// the number of "strands" is given by the number of chains + the number of patchy particles
	// since those are not linked to anything else
	*N_strands = _N_stars + _N_patchy;
	allocate_particles(particles);
}

void NathanStarInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(int i = 0; i < _N_patchy; i++) {
		NathanPatchyParticle *new_p = new NathanPatchyParticle();
		new_p->index = i;
		new_p->type = new_p->btype = PATCHY_PARTICLE;
		new_p->strand_id = i;

		particles[i] = new_p;
	}

	for(int i = _N_patchy; i < (int) particles.size(); i++) {
		NathanPolymerParticle *new_p = new NathanPolymerParticle();
		new_p->index = i;
		new_p->type = new_p->btype = POLYMER;
		new_p->strand_id = i;
		particles[i] = new_p;
	}
}

void NathanStarInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "T", &_T, 1);

	getInputNumber(&inp, "NATHAN_alpha", &_patch_alpha, 0);
	getInputNumber(&inp, "NATHAN_cosmax", &_patch_cosmax, 1);

	getInputInt(&inp, "NATHAN_interp_size", &_interp_size, 0);

	number size_ratio;
	getInputNumber(&inp, "NATHAN_size_ratio", &size_ratio, 1);
	_star_sigma_g = 1. / size_ratio;

	if(_star_sigma_g > 1.)
		_is_marzi = true;

	getInputInt(&inp, "NATHAN_f", &_star_f, 1);
	getInputNumber(&inp, "NATHAN_star_factor", &_star_factor, 0);

	getInputBool(&inp, "NATHAN_make_crystal", &_make_crystal, 0);
	if(_make_crystal) {
		getInputString(&inp, "NATHAN_crystal_type", _crystal_type, 1);
		_N_in_crystal = -1;
		getInputInt(&inp, "NATHAN_N_in_crystal", &_N_in_crystal, 0);
	}
}

void NathanStarInteraction::_setup_lambda_kappa() {
	if(_star_f < 2 || _star_f > 1e5)
		throw oxDNAException("NATHAN_f (%d) should be > 1 and < 10e5", _star_f);

	int fs[12] = { 2, 5, 10, 15, 18, 30, 40, 50, 65, 80, 100, 100000 };
	number lambdas[12] = { 0.46, 0.35, 0.30, 0.28, 0.27, 0.24, 0.24, 0.23, 0.23, 0.22, 0.22, 0.139 };
	number kappas[12] = { 0.58, 0.68, 0.74, 0.76, 0.77, 0.83, 0.85, 0.86, 0.87, 0.88, 0.89, 1. };

	int idx = 0;
	while(fs[idx] < _star_f)
		idx++;

	if(fs[idx] == _star_f) {
		_pi_lambda = lambdas[idx];
		_kappa_rs = kappas[idx];
	}
	else {
		number delta_x = (fs[idx] - fs[idx - 1]);
		number delta_f = (_star_f - fs[idx - 1]);
		number delta_l = lambdas[idx] - lambdas[idx - 1];
		number delta_k = kappas[idx] - kappas[idx - 1];

		_pi_lambda = lambdas[idx - 1] + delta_l / delta_x * delta_f;
		_kappa_rs = kappas[idx - 1] + delta_k / delta_x * delta_f;
	}
}

number NathanStarInteraction::_psi_1(number x, number smax) {
	return _xi * _exp_sqr_kappa_rs * (exp(-SQR(_kappa * x)) / x - exp(-SQR(_kappa * smax)) / smax) / _star_rs;
}

number NathanStarInteraction::_psi_2(number x, number smax) {
	return _xi * _exp_sqr_kappa_rs * (_sqrt_pi * (erf(_kappa * smax) - erf(_kappa * x)) / _kappa + x * exp(-SQR(_kappa * x)) - smax * exp(-SQR(_kappa * smax))) / _star_rs;
}

number NathanStarInteraction::_patchy_star_derivative(number r) {
	number sqr_r = r * r;
	number z = r - _patch_r;
	number smax = sqrt(z * (z + 1.));

	number derivative = -_T * _pi_lambda * _star_f3_2 * _patch_r / sqr_r;

	if(z < _star_rs)
		derivative *= (sqr_r - _sqr_patch_r) * (0.5 / SQR(z) - 0.5 / _sqr_star_rs + _psi_1(_star_rs, smax)) - log(z / _star_rs) + _psi_2(_star_rs, smax);
	else
		derivative *= (sqr_r - _sqr_patch_r) * _psi_1(z, smax) + _psi_2(z, smax);

	return derivative;
}

number NathanStarInteraction::_patchy_star_marzi_derivative(number r, gsl_spline *spl, gsl_interp_accel *acc, number bin) {
	number sqr_r = r * r;
	number z = r - _patch_r;
	number factor = -M_PI * _patch_r / sqr_r;

	if(spl == NULL) {
		spl = gsl_spline_alloc(gsl_interp_cspline, _interp_size);
		acc = gsl_interp_accel_alloc();

		double *ss = new double[_interp_size]();
		double *integrand = new double[_interp_size]();
		for(int i = 0; i < _interp_size; i++) {
			double s = (i + 0.5) * bin;
			double t = (z * (r + _patch_r) - SQR(s)) / s;
			ss[i] = s;
			integrand[i] = (sqr_r - _sqr_patch_r + SQR(s)) * (_pressure(s) - _pressure(s + t));
		}
		gsl_spline_init(spl, ss, integrand, _interp_size);

		delete[] ss;
		delete[] integrand;
	}

	number smax = sqrt(z * (z + 1.));
	number smax_used = (_interp_size - 0.5) * bin;
	if(smax > smax_used)
		smax = smax_used;
	if(z > smax)
		return 0.f;

	return factor * gsl_spline_eval_integ(spl, z, smax, acc);
}

number NathanStarInteraction::_pressure(number s) {
	number common_factor = _T * _pi_lambda * _star_f3_2 / M_PI;
	if(s < _star_rs)
		return common_factor / CUB(s);

	number sqr_s = SQR(s);
	number sqr_kappa = SQR(_kappa);
	return common_factor * (1. / sqr_s + 2 * sqr_kappa) * _xi * exp(-sqr_kappa * (sqr_s - _sqr_star_rs)) / _star_rs;
}

void NathanStarInteraction::_setup_interp() {
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
		if(_is_marzi)
			f_ys[i] = _patchy_star_marzi_derivative(r, marzi_spline, marzi_accel, bin);
		else
			f_ys[i] = _patchy_star_derivative(r);
	}
	gsl_spline_init(_spl_patchy_star, rs, f_ys, _interp_size);

	// print out the "real" derivative, its interpolation and their relative difference
	std::ofstream out_der("tabulated_patchy_star.dat");
	for(int i = 0; i < _interp_size - 1; i++) {
		r = rs[0] + 0.001 + bin * i;
		number real_der;
		if(_is_marzi)
			real_der = _patchy_star_marzi_derivative(r, marzi_spline, marzi_accel, bin);
		else
			real_der = _patchy_star_derivative(r);
		number interp_der = gsl_spline_eval(_spl_patchy_star, r, _acc_patchy_star);
		out_der << r << " " << real_der << " " << interp_der << " " << (real_der - interp_der) / real_der << std::endl;
	}
	out_der.close();

	double *u_ys = new double[_interp_size]();
	double u_shift = gsl_spline_eval_integ(_spl_patchy_star, rs[0], rs[_interp_size - 1], _acc_patchy_star);
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

void NathanStarInteraction::init() {
	_setup_lambda_kappa();

	_star_f1_2 = sqrt(_star_f);
	_star_f3_2 = CUB(_star_f1_2);

	_star_rg = 0.5 * _star_sigma_g;
	_star_rs = 2. * _star_rg / 3.;
	_star_sigma_s = 2. * _star_sigma_g / 3.;
	_sqr_star_rs = SQR(_star_rs);
	_sqr_star_sigma_s = SQR(_star_sigma_s);
	_sqr_star_sigma_g = SQR(_star_sigma_g);

	_patchy_star_rcut = _patch_r + _star_sigma_g;
	_sqr_patchy_star_rcut = SQR(_patchy_star_rcut);

	_star_rcut = 3. * _star_sigma_s;
	_sqr_star_rcut = SQR(_star_rcut);

	_kappa = _kappa_rs / _star_rs;
	_xi = 1. / (1. + 2. * SQR(_kappa_rs));
	_exp_sqr_kappa_rs = exp(SQR(_kappa_rs));
	_zeta = _sqrt_pi * _xi * erfc(_kappa_rs) * _exp_sqr_kappa_rs / _kappa_rs;
	_erfc_kappa_rs = erfc(_kappa_rs);

	_patch_cutoff = _patch_alpha * 1.5;
	_patchy_rcut = 1. + _patch_cutoff;
	_sqr_patchy_rcut = SQR(_patchy_rcut);
	this->_rcut = 1. + _patch_cutoff;
	if(_star_rcut > this->_rcut)
		this->_rcut = _star_rcut;
	this->_sqr_rcut = SQR(this->_rcut);
	_patch_angular_cutoff = _patch_cosmax * 1.5;

	_rep_E_cut = pow((number) this->_rcut, -_rep_power);

	_patch_pow_sigma = pow(_patch_cosmax, _patch_power);
	_patch_pow_alpha = pow(_patch_alpha, (number) 10.);

	_setup_interp();

	OX_LOG(Logger::LOG_INFO, "pi*lambda: %lf, kappa*rs: %lf, zeta: %lf, xi: %lf", _pi_lambda, _kappa_rs, _zeta, _xi);
	OX_LOG(Logger::LOG_INFO, "patchy rcut: %lf, patchy-star rcut: %lf, star-star rcut: %lf", _patchy_rcut, _patchy_star_rcut, _star_rcut);
}

void NathanStarInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

number NathanStarInteraction::_patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->type != PATCHY_PARTICLE && q->type != PATCHY_PARTICLE) {
		return 0.f;
	}
	number rnorm = _computed_r.norm();
	if(rnorm > _sqr_patchy_rcut) {
		return (number) 0.f;
	}

	// here everything is done as in Allen's paper
	number rmod = sqrt(rnorm);
	LR_vector r_versor = _computed_r / (-rmod);

	// repulsion
	number rep_part = 1. / pow(rnorm, _rep_power / 2);
	number energy = rep_part - _rep_E_cut;
	if(update_forces) {
		LR_vector force = r_versor * (_rep_power * rep_part / rmod);
		p->force += force;
		q->force -= force;
	}

	// attraction
	LR_vector p_axis = p->orientationT.v3;
	LR_vector q_axis = q->orientationT.v3;

	number cospr = -(p_axis * r_versor);
	if(cospr < 0.) {
		p_axis = -p_axis;
		cospr = -cospr;
	}
	number cosqr = q_axis * r_versor;
	if(cosqr < 0.) {
		q_axis = -q_axis;
		cosqr = -cosqr;
	}
	if(cospr < _patch_angular_cutoff || cosqr < _patch_angular_cutoff)
		return energy;

	number cospr_base = pow(cospr - 1., _patch_power - 1);
	// we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
	number cospr_part = cospr_base * (cospr - 1.);
	number p_mod = exp(-cospr_part / (2. * _patch_pow_sigma));

	number cosqr_base = pow(cosqr - 1., _patch_power - 1);
	number cosqr_part = cosqr_base * (cosqr - 1.);
	number q_mod = exp(-cosqr_part / (2. * _patch_pow_sigma));

	number sqr_surf_dist = SQR(rmod - 1.);
	number r8b10 = SQR(SQR(sqr_surf_dist)) / _patch_pow_alpha;
	number exp_part = -1.001 * exp(-(number) 0.5 * r8b10 * sqr_surf_dist);
	energy += exp_part * p_mod * q_mod;

	if(update_forces) {
		// radial part
		LR_vector tmp_force = r_versor * (p_mod * q_mod * 5. * (rmod - 1.) * exp_part * r8b10);

		// angular p part
		number der_p = exp_part * q_mod * (0.5 * _patch_power * p_mod * cospr_base / _patch_pow_sigma);
		LR_vector p_ortho = p_axis + cospr * r_versor;
		tmp_force -= p_ortho * (der_p / rmod);

		// angular q part
		number der_q = exp_part * p_mod * (-0.5 * _patch_power * q_mod * cosqr_base / _patch_pow_sigma);
		LR_vector q_ortho = q_axis - cosqr * r_versor;
		tmp_force -= q_ortho * (der_q / rmod);

		p->force += tmp_force;
		q->force -= tmp_force;

		p->torque += p->orientationT * (r_versor.cross(p_axis) * der_p);
		q->torque += q->orientationT * (r_versor.cross(q_axis) * der_q);
	}

	return energy;
}

number NathanStarInteraction::_patchy_star_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = (number) 0.f;
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_patchy_star_rcut) {
		return (number) 0.f;
	}

	number mod_r = sqrt(sqr_r);

	// this is just to avoid spitting out NaNs and it should occur very rarely, and only during equilibration
	if(mod_r < _spl_patchy_star->interp->xmin)
		energy = 1e6;
	// this can happen only in single precision
	else if(mod_r > _spl_patchy_star->interp->xmax)
		return (number) 0.f;
	else
		energy = gsl_spline_eval(_spl_patchy_star, mod_r, _acc_patchy_star);
	if(update_forces) {
		number force_mod = (mod_r < _spl_patchy_star->interp->xmin) ? 100 : -gsl_spline_eval_deriv(_spl_patchy_star, mod_r, _acc_patchy_star) / mod_r;
		p->force -= _computed_r * force_mod;
		q->force += _computed_r * force_mod;
	}
	return energy;
}

number NathanStarInteraction::_star_star_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = (number) 0.f;
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_star_rcut) {
		return (number) 0.f;
	}

	number mod_r = sqrt(sqr_r);

	number common_fact = _star_factor * 5. * _T * _star_f3_2 / 18.;
	number i_f = 1. / (1. + _star_f1_2 * 0.5);

	if(sqr_r < _sqr_star_sigma_s) {
		energy = -log(mod_r / _star_sigma_s) + i_f;

		if(update_forces) {
			// force over r
			number force_mod = common_fact / sqr_r;
			p->force -= _computed_r * force_mod;
			q->force += _computed_r * force_mod;
		}
	}
	else {
		number exp_factor = exp(-(mod_r - _star_sigma_s) * _star_f1_2 / (2. * _star_sigma_s));
		energy = i_f * (_star_sigma_s / mod_r) * exp_factor;

		if(update_forces) {
			// force over r
			number force_mod = common_fact * i_f * exp_factor * (_star_sigma_s / (sqr_r * mod_r) + _star_f1_2 / (2. * sqr_r));
			p->force -= _computed_r * force_mod;
			q->force += _computed_r * force_mod;
		}
	}

	return common_fact * energy;
}

number NathanStarInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number NathanStarInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.;
}

number NathanStarInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	int type = p->type + q->type;
	if(type == PATCHY_PATCHY) {
		return _patchy_interaction(p, q, false, update_forces);
	}
	else if(type == PATCHY_POLYMER) {
		return _patchy_star_interaction(p, q, false, update_forces);
	}
	else {
		return _star_star_interaction(p, q, false, update_forces);
	}
}

void NathanStarInteraction::generate_random_configuration(std::vector<BaseParticle*> &particles) {
	if(!_make_crystal) {
		BaseInteraction<NathanStarInteraction>::generate_random_configuration(particles);
		return;
	}
	// from here on it is implied that _make_crystal == true
	if(_N_in_crystal == -1)
		_N_in_crystal = _N_patchy;
	int N = particles.size();

	// check on the number of particles
	if(_N_in_crystal % 4 != 0)
		throw oxDNAException("The number of patchy particles should be a multiple of 4");
	int cxcycz = _N_in_crystal / 4;
	int cell_side = round(cbrt(cxcycz));
	if(CUB(cell_side) != cxcycz)
		throw oxDNAException("The number of patchy particles should be equal to 4*i^3, where i is an integer. For this specific case you may want to use %d", 4 * CUB(cell_side));

	OX_LOG(Logger::LOG_INFO, "Crystal cell_side: %d", cell_side);

	number gap = 0.05;
	number sigma_c = 1.;

	LR_vector tetrahedron_up[4];
	tetrahedron_up[0].x = tetrahedron_up[0].y = tetrahedron_up[0].z = 0;

	tetrahedron_up[1].x = sigma_c + gap;
	tetrahedron_up[1].y = tetrahedron_up[1].z = 0;

	tetrahedron_up[2].x = 0.5 * (sigma_c + gap);
	tetrahedron_up[2].y = sqrt(3.0) / 2.0 * (sigma_c + gap);
	tetrahedron_up[2].z = 0;

	tetrahedron_up[3].x = 0.5 * (sigma_c + gap);
	tetrahedron_up[3].y = 1. / 3. * sqrt(3.0) / 2.0 * (sigma_c + gap);
	tetrahedron_up[3].z = sqrt(2.0 / 3.0) * (sigma_c + gap);

	LR_matrix tetrahedron_orientT[4];
	for(int i = 0; i < 4; i++) {
		LR_vector axis(0, 0, 0);
		for(int j = 0; j < 4; j++) {
			if(j != i)
				axis += tetrahedron_up[j];
		}
		axis /= 3;
		axis -= tetrahedron_up[i];
		axis.normalize();
		tetrahedron_orientT[i] = Utils::get_random_rotation_matrix_from_angle(M_PI);
		tetrahedron_orientT[i].v3 = axis;
		tetrahedron_orientT[i].orthonormalize();
	}

	int cells[3] = { cell_side, cell_side, cell_side };

	// let's build the crystal
	if(_crystal_type == "hT") {
		LR_vector tetrahedron_down[4] = { tetrahedron_up[0], tetrahedron_up[1], tetrahedron_up[2], tetrahedron_up[3] };
		tetrahedron_down[2].y *= -1;
		tetrahedron_down[3].y *= -1;

		number xshift = (sigma_c + gap);
		number yshift = sqrt(3.0) * (sigma_c + gap);
		number zshift = 2.0 * sqrt(2.0 / 3.0) * (sigma_c + gap);

		LR_matrix tetrahedron_down_orientT[4];
		for(int i = 0; i < 4; i++) {
			LR_vector axis(0, 0, 0);
			for(int j = 0; j < 4; j++) {
				if(j != i)
					axis += tetrahedron_down[j];
			}
			axis /= 3;
			axis -= tetrahedron_down[i];
			axis.normalize();
			tetrahedron_down_orientT[i] = Utils::get_random_rotation_matrix_from_angle(M_PI);
			tetrahedron_down_orientT[i].v3 = axis;
			tetrahedron_down_orientT[i].orthonormalize();
		}

		int inserted = 0;
		for(int zc = 0; zc < cells[2]; zc++) {
			LR_vector *tetra = (zc % 2 == 0) ? tetrahedron_up : tetrahedron_down;
			LR_matrix *tetra_orientT = (zc % 2 == 0) ? tetrahedron_orientT : tetrahedron_down_orientT;
			number tot_zshift = zc * zshift;
			for(int yc = 0; yc < cells[1]; yc++) {
				number tot_yshift = yc * yshift;
				if(zc % 2 != 0)
					tot_yshift += yshift;
				for(int xc = 0; xc < cells[2]; xc++) {
					number tot_xshift = (yc % 2 + 2 * xc) * xshift;
					// loop over the tetrahedron
					for(int t = 0; t < 4; t++, inserted++) {
						BaseParticle *p = particles[inserted];
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
		number xshift = 2.0 * (sigma_c + gap);
		number yshift = 2.0 * sqrt(3.0) / 2.0 * (sigma_c + gap);
		number zshift = 2.0 * sqrt(2.0 / 3.0) * (sigma_c + gap);

		int inserted = 0;
		for(int zc = 0; zc < cells[2]; zc++) {
			number tot_zshift = 3 * zc * zshift;
			tot_zshift += (zc % 3) * zshift;
			tot_zshift = zc * zshift;

			for(int yc = 0; yc < cells[1]; yc++) {
				number tot_yshift = yc * yshift;
				if(zc % 3 == 1)
					tot_yshift -= 2. * yshift / 3.;
				else if(zc % 3 == 2)
					tot_yshift += 2. * yshift / 3.;
				for(int xc = 0; xc < cells[2]; xc++) {
					number tot_xshift = xc * xshift;
					if(yc % 2 == 1)
						tot_xshift += 0.5 * xshift;
					// loop over the tetrahedron
					for(int t = 0; t < 4; t++, inserted++) {
						BaseParticle *p = particles[inserted];
						p->pos.x = tetrahedron_up[t].x + tot_xshift;
						p->pos.y = tetrahedron_up[t].y + tot_yshift;
						p->pos.z = tetrahedron_up[t].z + tot_zshift;

						p->orientationT = tetrahedron_orientT[t];
						p->orientation = p->orientationT.get_transpose();
						p->set_positions();
						p->set_ext_potential(0, this->_box);
					}
				}
			}
		}
	}
	else
		throw oxDNAException("The crystal type be should be either 'hT' or 'cT' and not '%s'", _crystal_type.c_str());

	// now we add the polymers
	Cells c(particles, _box);
	c.init(_rcut);

	for(int i = _N_in_crystal; i < N; i++) {
		BaseParticle *p = particles[i];
		p->pos = LR_vector(0., 0., 0.);
	}

	c.global_update();

	LR_vector box_sides = this->_box->box_sides();
	for(int i = _N_in_crystal; i < N; i++) {
		BaseParticle *p = particles[i];

		bool inserted = false;
		do {
			p->pos = LR_vector(drand48() * box_sides.x, drand48() * box_sides.y, drand48() * box_sides.z);
			// random orientation
			p->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();

			p->set_positions();
			p->set_ext_potential(0, this->_box);
			c.single_update(p);

			inserted = true;

			// we take into account the bonded neighbours
			for(unsigned int n = 0; n < p->affected.size(); n++) {
				BaseParticle *p1 = p->affected[n].first;
				BaseParticle *p2 = p->affected[n].second;
				number e = 0.;
				if(p1->index <= p->index && p2->index <= p->index) {
					e = pair_interaction_bonded(p1, p2);
				}
				if(e > this->_energy_threshold)
					inserted = false;
			}

			// here we take into account the non-bonded interactions
			std::vector<BaseParticle*> neighs = c.get_complete_neigh_list(p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < p->index && this->generate_random_configuration_overlap(p, q))
					inserted = false;
			}

			// we take into account the external potential
			if(p->ext_potential > this->_energy_threshold) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && i % (N / 10) == 0)
			OX_LOG(Logger::LOG_INFO, "Inserted %d%% of the particles (%d/%d)", i*100/N, i, N);
	}
}

extern "C" NathanStarInteraction* make_NathanStarInteraction() {
	return new NathanStarInteraction();
}

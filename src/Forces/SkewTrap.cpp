/*
 * SkewTrap.cpp
 *
 *  Created on: 20/oct/2021
 *      Author: Jonah
 *  Updated on: 03/nov/2026
 *      Author: Erik
 */

#include "SkewTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"
#include <cmath>

namespace {

/**
 * Scaled complementary error function, erfcx(x) = exp(x^2) erfc(x), for x >= 0.
 *
 * Positive-argument branch of the double-precision implementation by njuffa (Shepherd & Laframboise,
 * Math. Comp. 36, 249 (1981)), https://stackoverflow.com/a/39777361, licensed under CC BY-SA 4.0.
 * Coefficients are written in decimal for C++14.
 * Maximum error is below 3 ulp. Negative arguments are not needed by SkewTrap and are not supported.
 */
double erfcx_nonneg(double x) {
	double a, d, e, m, p, q, r, s, t;

	a = std::fmax(x, 0.0 - x);

	// compute q = (a - 4) / (a + 4) accurately, mapping [0, inf) to [-1, 1]
	m = a - 4.0;
	p = a + 4.0;
	r = 1.0 / p;
	q = m * r;
	t = std::fma(q + 1.0, -4.0, a);
	e = std::fma(q, -a, t);
	q = std::fma(r, e, q);

	// approximate (1 + 2a) exp(a^2) erfc(a) as p(q) + 1 for q in [-1, 1]
	s = q * q;
	p = 8.9820305531190140e-10;
	t = 1.5764464777959401e-09;
	p = std::fma(p, s, -1.2155985739342269e-08);
	t = std::fma(t, s, -1.6386753783877791e-08);
	p = std::fma(p, s, 1.0585794011876720e-07);
	t = std::fma(t, s, 7.1190423171700940e-08);
	p = std::fma(p, s, -8.2040389712752056e-07);
	t = std::fma(t, s, 2.9796165315625938e-07);
	p = std::fma(p, s, 5.7059822144459833e-06);
	t = std::fma(t, s, -1.1225056665965572e-05);
	p = std::fma(p, s, -2.4397380523258482e-05);
	t = std::fma(t, s, 1.5062307184282616e-04);
	p = std::fma(p, s, -1.9925728768782324e-04);
	t = std::fma(t, s, -7.5777369791018515e-04);
	p = std::fma(p, s, 5.0319701025945277e-03);
	t = std::fma(t, s, -1.6197733983519948e-02);
	p = std::fma(p, s, 3.7167515521269866e-02);
	t = std::fma(t, s, -6.6330365820039094e-02);
	p = std::fma(p, s, 9.3732834999538536e-02);
	t = std::fma(t, s, -1.0103906603588378e-01);
	p = std::fma(p, s, 6.8097054254651804e-02);
	t = std::fma(t, s, 1.5379652102610957e-02);
	p = std::fma(p, q, t);
	p = std::fma(p, q, -1.3962111684056208e-01);
	p = std::fma(p, q, 2.3299511862555250e-01);

	// divide (1 + p) by (1 + 2a) to get exp(a^2) erfc(a)
	d = a + 0.5;
	r = 1.0 / d;
	r = r * 0.5;
	q = std::fma(p, r, r);
	t = q + q;
	e = (p - q) + std::fma(t, -a, 1.0);
	r = std::fma(e, r, q);

	if(a > 1.7976931348623157e+308) r = 0.0; // infinite argument

	return r;
}

}

SkewTrap::SkewTrap() :
				BaseForce() {
	_ref_id = -2;
	_particle = -2;
	_p_ptr = NULL;
	_r0 = -1.;
	_a = 0.;
	_s = 0.;
	PBC = false;
	_inv_s2 = 0.;
	_w_coeff = 0.;
	_g_coeff = 0.;
	_V_const = 0.;
}

std::tuple<std::vector<int>, std::string> SkewTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particle_string;
	std::string ref_particle_string;
	getInputString(&inp, "particle", particle_string, 1);
	getInputString(&inp, "ref_particle", ref_particle_string, 1);
	_particle = Utils::get_single_particle_from_string(CONFIG_INFO->particles(), particle_string, "SkewTrap particle");
	_ref_id = Utils::get_single_particle_from_string(CONFIG_INFO->particles(), ref_particle_string, "SkewTrap ref_particle");
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stdev", &_s, 1);
	getInputNumber(&inp, "shape", &_a, 1);
	if(_s <= 0.) {
		throw oxDNAException("SkewTrap: stdev must be positive (got %g)", _s);
	}
	if(_a == 0.) {
		throw oxDNAException("SkewTrap cannot run with a shape parameter of 0, please use mutual_trap");
	}
	getInputBool(&inp, "PBC", &PBC, 0);
	_rate = 0.; //default rate is 0
	getInputNumber(&inp, "rate", &_rate, 0);

	_inv_s2 = 1. / SQR(_s);
	_w_coeff = -_a / (_s * sqrt(2.));
	_g_coeff = _a * _s * sqrt(2. / PI);
	_V_const = log(sqrt(2. * PI) * _s);

	_p_ptr = CONFIG_INFO->particles()[_ref_id];

	std::string description = Utils::sformat("SkewTrap (stdev=%g, shape=%g, rate=%g, r0=%g, ref_particle=%d, PBC=%d)", _s, _a, _rate, _r0, _ref_id, PBC);

	return std::make_tuple(std::vector<int>{_particle}, description);
}

LR_vector SkewTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC) {
		return CONFIG_INFO->box->min_image(u, v);
	}
	else {
		return v - u;
	}
}

number SkewTrap::potential_of_x(number x) const {
	number w = _w_coeff * x;
	// ln erfc(w), using erfc(w) = erfcx(w) exp(-w^2) on the steep side, where erfc(w) is tiny
	number log_erfc = (w > 0.) ? log(erfcx_nonneg(w)) - SQR(w) : log(erfc(w));
	return 0.5 * SQR(x) * _inv_s2 + _V_const - log_erfc;
}

number SkewTrap::force_mag_of_x(number x) const {
	number w = _w_coeff * x;
	// g = exp(-w^2) / erfc(w)
	number g = (w > 0.) ? 1. / erfcx_nonneg(w) : exp(-SQR(w)) / erfc(w);
	return (-x + _g_coeff * g) * _inv_s2;
}

LR_vector SkewTrap::force(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
	number dr_abs = dr.module();
	number x = dr_abs - (_r0 + (_rate * step));
	return -(dr / dr_abs) * force_mag_of_x(x);
}

number SkewTrap::potential(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
	number x = dr.module() - (_r0 + (_rate * step));
	return potential_of_x(x);
}

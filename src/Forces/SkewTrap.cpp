/*
 * SkewTrap.cpp
 *
 *  Created on: 20/oct/2021
 *      Author: Jonah
 */

#include "SkewTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"
#include <cmath>

SkewTrap::SkewTrap() :
				BaseForce() {
	_ref_id = -2;
	_particle = -2;
	_p_ptr = NULL;
	_r0 = -1.;
	PBC = false;
	_val1 = 0.f;
	_val2 = 0.f;
	_val3 = 0.f;
	_val4 = 0.f;
	_val5 = 0.f;
	_val6 = 0.f;
	_a = 0.f;
	_s = 0.f;
}

std::tuple<std::vector<int>, std::string> SkewTrap::init(input_file &inp) {
	getInputInt(&inp, "particle", &_particle, 1);
	getInputInt(&inp, "ref_particle", &_ref_id, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stdev", &_s, 1);
	getInputNumber(&inp, "shape", &_a, 1);
    if(_a == 0.f) throw oxDNAException("Skew Force Cannot Run with a shape parameter of 0, please use mutual trap");
	getInputBool(&inp, "PBC", &PBC, 0);
	_rate = 0.f; //default rate is 0
	getInputNumber(&inp, "rate", &_rate, 0);

	// set precomputed values
	_val1 = 1.f / (2.f*pow(_s, 2)); // 1/(2s^2)
	_val2 = sqrt(2.f*PI) * _s; // Sqrt(2*Pi) * s
	_val3 = -1.f*pow(_a, 2)/ (2*pow(_s, 2)); // - a^2/(2s^2)
	_val4 = _a/(_s*sqrt(2.f)); // a/ (s*sqrt(2))
    _val5 = 1.f / (pow(_s, 2)); // 1/(s^2)
    _val6 = sqrt(2/PI) * _s; //sqrt(2/pi) *s

    _ddx = 0.f;
    _slope = 0.f;
    _intercept = 0.f;

    _find_discontinuity(); // fills the 3 values above for our fit function to avoid discontinuities in force

	int N = CONFIG_INFO->particles().size();
	if(_ref_id < 0 || _ref_id >= N) {
		throw oxDNAException("Invalid reference particle %d for Skew Trap", _ref_id);
	}
	_p_ptr = CONFIG_INFO->particles()[_ref_id];

	if(_particle >= N || N < -1) {
		throw oxDNAException("Trying to add a SkewTrap on non-existent particle %d. Aborting", _particle);
	}
	if(_particle == -1) {
		throw oxDNAException("Cannot apply SkewTrap to all particles. Aborting");
	}

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

LR_vector SkewTrap::value(llint step, LR_vector &pos) {
    // Calculates: (- x + (a e^((-1/2s^2) *a^2 * x^2) * Sqrt(2/Pi)) / (1 + Erf[(a*x)/Sqrt(2)])]) /  (s^2)
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
    number dx = dr.module() - (_r0 + (_rate * step));
    number f_mag;
    if((dx > 0) == (_ddx > 0) && (abs(dx) >= abs(_ddx))){ // same sign and past discontinuous point
        // Hep F (don't worry it's not contagious)
        f_mag = -8.f * _slope * pow(dx, 7) - _intercept;
        //printf("octa\n");
    } else {
        number numerator = _a*exp(pow(dx, 2) * _val3) * _val6; // (a e^((-a^2/(2s^2))* x^2) * Sqrt(2/Pi)*s)
        number denominator = 1 + erf(dx*_val4); // 1 +Erf[ (x) * a / (s*Sqrt(2)) ]
        f_mag = (-dx + numerator/denominator)*_val5;
        //printf("skewforce\n");
    }
    //debug
    //printf("dx %.4f", dx);
    //printf("force: %.4f \n", f_mag);
    return -(dr / dr.module()) * f_mag;
}

number SkewTrap::potential(llint step, LR_vector &pos) {
    // Calculates: Log[(e^(x^2/(2s^2))*Sqrt(2Pi)*s)  /  (1 + Erf[(a*x)/(s*Sqrt(2)]))]
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
    number dx = dr.module() - (_r0 + (_rate * step));

    number pot;

    if((dx > 0) == (_ddx > 0) && (abs(dx) >= abs(_ddx))){ // same sign and past discontinuous point
        // Oct E
        pot = _slope*pow(dx, 8) + _intercept * dx;
    } else {
        //normal skew value
        number numerator = exp(pow(dx, 2) * _val1) * _val2; // e^(x^2 *(1/(2s^2))*Sqrt(2Pi)) *s)
        number denominator = 1 + erf(dx*_val4); // 1 +Erf[ (x) * a / (s*Sqrt(2)) ]
        pot = log(numerator/denominator);
    }

	return pot;
}

number SkewTrap::_solution_term(number dx) {
    // Calculates: Log[(e^(x^2/(2s^2))*Sqrt(2Pi)*s)  /  (1 + Erf[(a*x)/(s*Sqrt(2)]))]
    number numerator = exp(pow(dx, 2) * _val1) * _s; // e^(x^2 *(1/(2s^2))*Sqrt(2Pi)) *s)
    number denominator = 1 + erf(dx*_val4); // 1 +Erf[ (x) * a / (s*Sqrt(2)) ]
    number pot = log(numerator/denominator);
    return pot;
}

// helper function for discontinuity calculation
number SkewTrap::_force_mag(number dx) {
    // Calculates: (- x + (a e^((-1/2s^2) *a^2 * x^2) * Sqrt(2/Pi)) / (1 + Erf[(a*x)/Sqrt(2)])]) /  (s^2)
    number numerator = _a*exp(pow(dx, 2) * _val3) * _val6; // (a e^((-a^2/(2s^2))* x^2) * Sqrt(2/Pi)*s)
    number denominator = 1 + erf(dx*_val4); // 1 +Erf[ (x) * a / (s*Sqrt(2)) ]
    number f = -(-dx + numerator/denominator)*_val5;;
    return f;
}

void SkewTrap::_find_discontinuity()
{
    number x = 0.02;
    bool discontinuity_found = false;
    number oldforcepos = _force_mag(x);
    number oldforceneg = _force_mag(-x);
    number forcepos = _force_mag(x+0.002f), forceneg = _force_mag(-x-0.002f);
    number posdir = (forcepos - oldforcepos)/abs(forcepos - oldforcepos);
    number negdir = (forceneg - oldforceneg)/abs(forceneg - oldforceneg);


    // discontinuity can be found by finding where the slope of the function changes sign
    // discontinuity can occur on either side (+ or -) depending on shape parameter
    // added buffer space of 0.004 from discontinuity point
    while(!discontinuity_found){
        x += 0.002;
        forcepos = _force_mag(x);
        forceneg = _force_mag(-x);

        // check if sign of the slope of the force function changes
        if((forcepos - oldforcepos <0) == (posdir<0)){
            oldforcepos = forcepos;
        }else{
            discontinuity_found = true;
            _ddx = x - 0.004;
            _calc_slope();
            _calc_intercept();
        }

        // check if sign of the slope of the force function changes
        if((forceneg - oldforceneg <0) == (negdir<0)){
            oldforceneg = forceneg;
        }else{
            discontinuity_found = true;
            _ddx = -x + 0.004;
            _calc_slope();
            _calc_intercept();
        }
    }

    // Debug
    // _ddx = 0.4f;
    // _calc_slope();
    // _calc_intercept();
     printf("Discontinuity x %.4f slope %.4f intercept %.4f \n", _ddx, _slope, _intercept);
     printf("Skew Params s %.4f a %.4f \n", _s, _a);
};

void SkewTrap::_calc_slope(){
    // oct E hep F
    number term1 = _ddx*((_a*exp(- pow(_a, 2)* pow(_ddx,2)/ (2.f*pow(_s, 2)))*_val6)/(1+erf(_ddx*_val4)) -_ddx) / SQR(_s);
    number numerator = exp(pow(_ddx, 2) * _val1) * _val2; // e^(x^2 *(1/(2s^2))*Sqrt(2Pi)) *s)
    number denominator = 1 + erf(_ddx*_val4); // 1 +Erf[ (x) * a / (s*Sqrt(2)) ]
    number term2 = log(numerator/denominator);
    _slope = -(term1 + term2)/(7.f * pow(_ddx, 8));
}

void SkewTrap::_calc_intercept(){
    // oct E hep F
    number term1 = _ddx/pow(_s, 2); // 2 / s^2
    number term2 = (_a*exp(- pow(_a, 2)* pow(_ddx,2)/ (2.f*pow(_s, 2))) * sqrt(2.f/PI)) /
                   (_s*(1+erf(_ddx*_val4)));
    number term3 = 4.f*log(2*PI)/_ddx;
    number term4 = (8*_solution_term(_ddx))/_ddx;

    _intercept = (-term1 + term2 + term3 + term4)/7.f;
}

/* The code in the below section is kept in case we change function that is fit to go past the non continuity in the forces
 *
 *
 *
 * #### Potential
 * Quadratic E
   pot = _slope/2 * pow(dx, 2) + _intercept * dx;

   Quart E
   pot = _slope * pow(dx, 4) + _intercept * dx;
 *
 * #### Force
 *  Lin F
    f_mag = -_ddx*_slope - _intercept;

    Cub F
    f_mag = -4.f*pow(_ddx, 3)*_slope - _intercept;

 * #### Calc Intercept
 *  ##quadratic E Linear F
    number term1 = _ddx/pow(_s, 2); // 2 / s^2
    number term2 = (_a*exp(- pow(_a, 2)* pow(_ddx,2)/ (2.f*pow(_s, 2))) * sqrt(2.f/PI)) /
                   (_s*(1+erf(_ddx*_val4)));
    number term3 = log(2*PI)/_ddx;
    number term4 = (2*_solution_term(_ddx))/_ddx;
    _intercept = -term1 + term2 + term3 + term4;

    ##quart E cub F
    number term1 = _ddx/pow(_s, 2); // x / s^2
    number term2 = (_a*exp(- pow(_a, 2)* pow(_ddx,2)/ (2.f*pow(_s, 2))) * sqrt(2.f/PI)) /
                   (_s*(1+erf(_ddx*_val4)));
    number term3 = 2.f*log(2*PI)/_ddx;
    number term4 = (4*_solution_term(_ddx))/_ddx;

    _intercept = (-term1 + term2 + term3 + term4)/3.f;
 *
 * ### Calc Slope
 *  ##quadratic E Linear F
    number term1 = 2.f/pow(_s, 2); // 2 / s^2
    number term2 = (2.f*_a*exp(- pow(_a, 2)* pow(_ddx,2)/ (2.f*pow(_s, 2))) * sqrt(2.f/PI)) /
            (_s*_ddx*(1+erf(_ddx*_val4)));
    number term3 = (log(2*PI) + 2* _solution_term(_ddx))/SQR(_ddx);
    _slope = term1 - term2 - term3;

    ##quart E cub F
    number term1 = _ddx*((_a*exp(- pow(_a, 2)* pow(_ddx,2)/ (2.f*pow(_s, 2)))*_val6)/(1+erf(_ddx*_val4)) -_ddx) / SQR(_s);
    number numerator = exp(pow(_ddx, 2) * _val1) * _val2; // e^(x^2 *(1/(2s^2))*Sqrt(2Pi)) *s)
    number denominator = 1 + erf(_ddx*_val4); // 1 +Erf[ (x) * a / (s*Sqrt(2)) ]
    number term2 = log(numerator/denominator);
    _slope = -(term1 + term2)/(3.f * pow(_ddx, 4));
 */




















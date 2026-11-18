/*
 * MaxBackboneForce.cpp
 */

#include "MaxBackboneForce.h"

#include "ConfigInfo.h"
#include "Logger.h"
#include "oxDNAException.h"
#include "parse_input/parse_input.h"

#include <cmath>

MaxBackboneForce::MaxBackboneForce(FarForce far_force) :
				_far_force_mode(far_force) {

}

void MaxBackboneForce::get_settings(input_file &inp) {
	if(getInputNumber(&inp, "max_backbone_force", &_fmax, 0) != KEY_FOUND) {
		return;
	}

	_enabled = true;
	if(_fmax < 0.f) {
		throw oxDNAException("Cowardly refusing to run with a negative max_backbone_force");
	}
	_fmax_initial = _fmax;

	getInputNumber(&inp, "max_backbone_rate", &_rate, 0);
	if(_rate < (number) 0.) {
		throw oxDNAException("Cowardly refusing to run with a negative max_backbone_rate");
	}

	if(_far_force_mode == FarForce::INPUT) {
		if(getInputNumber(&inp, "max_backbone_force_far", &_finf, 0) != KEY_FOUND) {
			_finf = 0.04f; // roughly 2pN, very weak but still there
		}
		if(_finf < 0.f) {
			throw oxDNAException("Cowardly refusing to run with a negative max_backbone_force_far");
		}
	}
	else {
		number tmp;
		if(getInputNumber(&inp, "max_backbone_force_far", &tmp, 0) == KEY_FOUND) {
			OX_LOG(Logger::LOG_WARNING, "max_backbone_force_far is ignored by this interaction: the far force is set by max_backbone_force and the FENE parameters");
		}
	}
}

void MaxBackboneForce::set_fene(number eps, number delta2) {
	_eps = eps;
	_delta2 = delta2;
	_fene_set = true;
	if(_enabled) {
		_compute();
	}
}

void MaxBackboneForce::log_settings(const char *interaction_name) const {
	if(!_enabled) {
		return;
	}
	OX_LOG(Logger::LOG_INFO, "%s: using a maximum backbone force of %g (the corresponding mbf_xmax is %g) and a far value of %g", interaction_name, _fmax, _xmax, _finf);
	if(_rate > (number) 0.) {
		OX_LOG(Logger::LOG_INFO, "%s: the maximum backbone force grows by %g per step (max_backbone_rate)", interaction_name, _rate);
	}
}

bool MaxBackboneForce::update() {
	if(!_enabled || _rate == (number) 0.) {
		return false;
	}
	llint step = CONFIG_INFO->curr_step;
	if(step == _step) {
		return false;
	}
	_step = step;
	_fmax = _fmax_initial + _rate * step;
	if(_fene_set) {
		_compute();
	}
	return true;
}

void MaxBackboneForce::_compute() {
	_xmax = xmax(_fmax, _eps, _delta2);
	if(_far_force_mode == FarForce::CONTINUOUS) {
		// the far force that makes the capped potential continuous with the FENE one without a constant
		number fene_xmax = -(_eps / 2.f) * log(1.f - SQR(_xmax) / _delta2);
		number log_coefficient = (fene_xmax - _fmax * _xmax) / (log(_xmax) - 1.0);
		_finf = _fmax - log_coefficient / _xmax;
	}
}

number MaxBackboneForce::offset() const {
	number fene_xmax = -(_eps / 2.f) * log(1.f - SQR(_xmax) / _delta2);
	return fene_xmax - (_fmax - _finf) * _xmax * log(_xmax) - _finf * _xmax;
}

number MaxBackboneForce::xmax(number fmax, number eps, number delta2) {
	return (-eps + sqrt(eps * eps + 4.f * fmax * fmax * delta2)) / (2.f * fmax);
}

number MaxBackboneForce::energy(number x, number fmax, number finf, number xmax, number eps, number delta2) {
	number fene_xmax = -(eps / 2.f) * log(1.f - SQR(xmax) / delta2);
	number long_xmax = (fmax - finf) * xmax * log(xmax) + finf * xmax;
	return (fmax - finf) * xmax * log(fabs(x)) + finf * fabs(x) - long_xmax + fene_xmax;
}

number MaxBackboneForce::force(number x, number fmax, number finf, number xmax) {
	return -copysign(1.f, x) * ((fmax - finf) * xmax / fabs(x) + finf);
}

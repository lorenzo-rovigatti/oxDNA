/*
 * GenericCentralForce.cpp
 *
 *  Created on: 22/jun/2017
 *      Author: Lorenzo
 */

#include "GenericCentralForce.h"

#include <cfloat>
#include <fstream>
#include <map>

using namespace std;

LookupTable::LookupTable() {
	_N = -1;
	_xupp = -1.;
	_xlow = -1.;
	_inv_sqr_delta = -1.;
	_delta = -1.;
}

LookupTable::~LookupTable() {

}

number LookupTable::_linear_interpolation(number x, vector<number> &x_data, vector<number> &fx_data) {
	int ind = -1;
	for(uint32_t i = 0; i < x_data.size() && ind == -1; i++) {
		if(x_data[i] > x) ind = i;
	}

	number val;
	if(ind == -1) {
		int last = fx_data.size() - 1;
		number slope = (fx_data[last] - fx_data[last - 1]) / (x_data[last] - x_data[last - 1]);
		val = fx_data[last] + slope * (x - x_data[last]);
	}
	else {
		number slope = (fx_data[ind] - fx_data[ind - 1]) / (x_data[ind] - x_data[ind - 1]);
		val = fx_data[ind - 1] + slope * (x - x_data[ind - 1]);
	}

	return val;
}

void LookupTable::init(std::string filename, int N) {
	_N = N;

	ifstream input(filename.c_str());
	if(!input.good()) throw oxDNAException("GenericCentralForce: potential_file '%s' not found", filename.c_str());

	vector<number> v_x, v_fx, v_dfx;
	while(input.good()) {
		number x, fx, dfx;
		input >> x >> fx >> dfx;

		if(!input.fail()) {
			if(v_x.size() > 0 && x <= v_x.back()) throw oxDNAException("GenericCentralForce: the x axis in the potential file should be sorted in an ascending fashion (%lf > %lf)", v_x.back(), x);

			v_x.push_back(x);
			v_fx.push_back(fx);
			v_dfx.push_back(dfx);
		}
	}
	input.close();

	_xlow = v_x.front();
	_xupp = v_x.back();

	_delta = (_xupp - _xlow) / (number) _N;
	_inv_sqr_delta = 1. / SQR(_delta);
	OX_LOG(Logger::LOG_INFO, "GenericCentralForce: detected x_low = %lf, x_upp = %lf", _xlow, _xupp);

	_A.reserve(_N + 1);
	_B.reserve(_N + 1);
	_C.reserve(_N + 1);
	_D.reserve(_N + 1);
	for(int i = 0; i <= _N; i++) {
		number x = _xlow + i * _delta;

		number fx0 = _linear_interpolation(x, v_x, v_fx);
		number fx1 = _linear_interpolation(x + _delta, v_x, v_fx);
		number derx0 = _linear_interpolation(x, v_x, v_dfx);
		number derx1 = _linear_interpolation(x + _delta, v_x, v_dfx);

		_A[i] = fx0;
		_B[i] = derx0;
		_D[i] = (2 * (fx0 - fx1) + (derx0 + derx1) * _delta) / _delta;
		_C[i] = (fx1 - fx0 + (-derx0 - _D[i]) * _delta);
	}
}

inline number LookupTable::query_function(number x) {
	if(x <= _xlow) return _A[0];
	if(x >= _xupp) return _A[_N];
	int i = (int) ((x - _xlow) / _delta);
	number dx = x - _xlow - _delta * i;
	return (_A[i] + dx * (_B[i] + dx * (_C[i] + dx * _D[i]) * _inv_sqr_delta));
}

inline number LookupTable::query_derivative(number x) {
	if(x < _xlow) return _B[0];
	if(x >= _xupp) return _B[_N];
	int i = (int) ((x - _xlow) / _delta);
	number dx = x - _xlow - _delta * i;
	return (_B[i] + (2 * dx * _C[i] + 3 * SQR(dx) * _D[i]) * _inv_sqr_delta);
}


GenericCentralForce::GenericCentralForce() :
				BaseForce() {
	_type = -1;
	_table_N = -1;
	_E_shift = 0.;

	inner_cut_off = inner_cut_off_sqr = 0.;
	outer_cut_off = outer_cut_off_sqr = 0.;
	_supported_types[string("gravity")] = GRAVITY;
	_supported_types[string("interpolated")] = INTERPOLATED;
}

GenericCentralForce::~GenericCentralForce() {

}

std::tuple<std::vector<int>, std::string> GenericCentralForce::init(input_file &inp) {
	BaseForce::init(inp);

	string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "E_shift", &_E_shift, 0);

	string strdir;
	getInputString(&inp, "center", strdir, 1);
	vector<string> spl = Utils::split(strdir, ',');
	if(spl.size() != 3) throw oxDNAException("Could not parse 'center' in external_forces_file. Dying badly");

	center.x = atof(spl[0].c_str());
	center.y = atof(spl[1].c_str());
	center.z = atof(spl[2].c_str());

	string type_str;
	getInputString(&inp, "force_type", type_str, 1);
	map<string, int>::iterator res = _supported_types.find(type_str);
	if(res == _supported_types.end()) throw oxDNAException("GenericCentralForce: unsupported type '%s'", type_str.c_str());

	_type = res->second;
	switch(_type) {
	case GRAVITY:
		getInputNumber(&inp, "F0", &_F0, 1);
		getInputNumber(&inp, "inner_cut_off", &inner_cut_off, 0);
		getInputNumber(&inp, "outer_cut_off", &outer_cut_off, 0);
		break;
	case INTERPOLATED:
		getInputString(&inp, "potential_file", _table_filename, 1);
		getInputInt(&inp, "interpolated_N", &_table_N, 1);
		break;
	default:
		break;
	}

	inner_cut_off_sqr = SQR(inner_cut_off);
	outer_cut_off_sqr = SQR(outer_cut_off);

	switch(_type) {
	case INTERPOLATED:
		_table.init(_table_filename, _table_N);
		break;
	default:
		break;
	}

	std::string description = Utils::sformat("GenericCentralForce (center=%g,%g,%g)", center.x, center.y, center.z);
	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "ConstantRateForce");

	return std::make_tuple(particle_ids, description);
}

LR_vector GenericCentralForce::value(llint step, LR_vector &pos) {
	LR_vector dir = center - pos;
	number dist_sqr = dir.norm();
	dir.normalize();

	switch(_type) {
	case GRAVITY:
		if(dist_sqr < inner_cut_off_sqr) {
			return LR_vector(0., 0., 0.);
		}
		if(outer_cut_off > 0. && dist_sqr > outer_cut_off_sqr) {
			return LR_vector(0., 0., 0.);
		}

		return _F0 * dir;
	case INTERPOLATED:
		return _table.query_derivative(sqrt(dist_sqr)) * dir;
	}

	// we should never reach this point, but we return something to remove warnings
	return LR_vector(0., 0., 0.);
}

number GenericCentralForce::potential(llint step, LR_vector &pos) {
	number dist = (center - pos).module();
	number energy = 0.f;

	switch(_type) {
	case GRAVITY:
		if(dist < inner_cut_off) {
			energy = 0.;
		}
		else if(outer_cut_off > 0. && dist > outer_cut_off) {
			energy = 0.;
		}
		else {
			energy = _F0 * dist - _F0 * inner_cut_off + _E_shift;
		}
		break;
	case INTERPOLATED:
		energy = _table.query_function(dist);
		break;
	}
	return energy;
}


/*
 * PolydisperseLTInteraction.cpp
 *
 *  Created on: 20 Jun 2018
 *      Author: lorenzo
 */

#include "PolydisperseLTInteraction.h"
#include "Particles/CustomParticle.h"

PolydisperseLTInteraction::PolydisperseLTInteraction() {
	_lt_points = 100;
	_is_polydisperse = true;

	ADD_INTERACTION_TO_MAP(BONDED, pair_interaction_bonded);
	ADD_INTERACTION_TO_MAP(NONBONDED, pair_interaction_nonbonded);
}

PolydisperseLTInteraction::~PolydisperseLTInteraction() {

}

number PolydisperseLTInteraction::_linear_interpolation(number x, number *x_data, number *fx_data, int points) {
	int ind = -1;
	for(int i = 0; i < points && ind == -1; i++)
		if(x_data[i] > x) ind = i;

	number val;
	if(ind == -1) {
		int last = points - 1;
		number slope = (fx_data[last] - fx_data[last - 1]) / (x_data[last] - x_data[last - 1]);
		val = fx_data[last] + slope * (x - x_data[last]);
	}
	else {
		number slope = (fx_data[ind] - fx_data[ind - 1]) / (x_data[ind] - x_data[ind - 1]);
		val = fx_data[ind - 1] + slope * (x - x_data[ind - 1]);
	}

	return val;
}

number PolydisperseLTInteraction::_fx(number x, void *par) {
	base_function *data = (base_function *) par;
	return _linear_interpolation(x, data->x, data->fx, data->points);
}

number PolydisperseLTInteraction::_dfx(number x, void *par) {
	base_function *data = (base_function *) par;
	return _linear_interpolation(x, data->x, data->dfx, data->points);
}

void PolydisperseLTInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputString(&inp, "custom_lt_file", _lt_filename, 1);
	getInputInt(&inp, "custom_points", &_lt_points, 0);
	getInputBool(&inp, "custom_is_polydisperse", &_is_polydisperse, 0);
	getInputNumber(&inp, "custom_rcut", &_rcut_base, 1);
}

void PolydisperseLTInteraction::init() {
	std::ifstream lt_file(_lt_filename, std::ios::in);
	if(!lt_file.good()) throw oxDNAException("Can't read lookup file '%s'. Aborting", _lt_filename);

	std::string line;
	int n_lines = 0;
	bool stop = false;
	while(!stop) {
		getline(lt_file, line);
		std::vector<std::string> spl = Utils::split(line);
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else n_lines++;
	}

	if(n_lines == 0) throw oxDNAException("Zero lines found in the lookup file. The only parsed lines are those with three columns (x, U, F)");

	base_function data;
	data.x = new number[n_lines];
	data.fx = new number[n_lines];
	data.dfx = new number[n_lines];
	data.points = n_lines;

	lt_file.clear();
	lt_file.seekg(0, std::ios::beg);
	int i = 0;
	stop = false;
	while(!stop) {
		getline(lt_file, line);
		std::vector<std::string> spl = Utils::split(line);
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else {
			data.x[i] = atof(spl[0].c_str());
			data.fx[i] = atof(spl[1].c_str());
			data.dfx[i] = atof(spl[2].c_str());

			if(i > 0 && data.x[i] <= data.x[i - 1]) throw oxDNAException("The x values of the lookup table should be monotonically increasing (found x[%d] = %f <= %f = x[%d])", i, data.x[i], i - 1, data.x[i - 1]);
			i++;
		}
	}
	lt_file.close();

	number lowlimit = data.x[0];
	number uplimit = data.x[i - 1];

	_lookup_table.build([this](number x, void *args) { return this->_fx(x, args); }, [this](number x, void *args) { return this->_dfx(x, args); }, (void *) (&data), _lt_points, lowlimit, uplimit);

	delete[] data.x;
	delete[] data.fx;
	delete[] data.dfx;

	_rcut = _rcut_base;
	_Ecut = _lookup_table.query(_rcut_base);
}

void PolydisperseLTInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(int i = 0; i < (int) particles.size(); i++)
		particles[i] = new CustomParticle();
}

void PolydisperseLTInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	*N_strands = N;

	std::ifstream topology(_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[512];
	topology.getline(line, 512);

	_sigmas.resize(N, 0.);

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = particles[i]->strand_id = i;
		if(_is_polydisperse) {
			topology >> _sigmas[i];
			if(_sigmas[i] <= (number) 0.) throw oxDNAException("The %d-th particle has either no associated diameter or the diameter is <= 0", i);
			if(_sigmas[i] > _rcut) _rcut = _rcut_base * _sigmas[i];
		}
	}
	topology.close();

	_sqr_rcut = SQR(_rcut);
	OX_LOG(Logger::LOG_INFO, "custom: rcut = %lf, Ecut = %lf", _rcut, _Ecut);
}

number PolydisperseLTInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number PolydisperseLTInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return 0.;
}

number PolydisperseLTInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) return 0.;

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	number sigma = (_is_polydisperse) ? (_sigmas[p->index] + _sigmas[q->index]) / 2. : 1;
	_computed_r /= sigma;

	number dist = _computed_r.module();
	if(dist > _rcut_base) return 0.;

	number energy = _lookup_table.query(dist) - _Ecut;

	if(dist < _lookup_table.x_low()) {
		fprintf(stderr, "Exceeded the lower bound (%lf < %lf)\n", dist, _lookup_table.x_low());
	}

	if(update_forces) {
		number force_mod = -_lookup_table.query_derivative(dist) / sigma;
		LR_vector force = _computed_r * (force_mod / dist);
		p->force -= force;
		q->force += force;
	}

	return energy;
}

void PolydisperseLTInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

extern "C" PolydisperseLTInteraction *make_PolydisperseLTInteraction() {
	return new PolydisperseLTInteraction();
}

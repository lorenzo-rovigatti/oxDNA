/*
 * CustomInteraction.cpp
 *
 *  Created on: 6 Mar 2014
 *      Author: lorenzo
 */

#include "CustomInteraction.h"
#include "../Particles/CustomParticle.h"

template<typename number>
CustomInteraction<number>::CustomInteraction() {
	_bonded_points = 100;

	this->_int_map[BONDED] = &CustomInteraction<number>::pair_interaction_bonded;
	this->_int_map[NONBONDED] = &CustomInteraction<number>::pair_interaction_nonbonded;
}

template<typename number>
CustomInteraction<number>::~CustomInteraction() {

}

template<typename number>
number CustomInteraction<number>::_linear_interpolation(number x, number *x_data, number *fx_data, int points) {
	int ind = -1;
	for(int i = 0; i < points && ind == -1; i++) if(x_data[i] > x) ind = i;

	number val;
	if(ind == -1) {
		int last = points - 1;
		number slope = (fx_data[last] - fx_data[last-1]) / (x_data[last] - x_data[last-1]);
		val = fx_data[last] + slope * (x - x_data[last]);
	}
	else {
		number slope = (fx_data[ind] - fx_data[ind-1]) / (x_data[ind] - x_data[ind-1]);
		val = fx_data[ind-1] + slope * (x - x_data[ind-1]);
	}

	return val;
}

template<typename number>
number CustomInteraction<number>::_fx(number x, void *par) {
	base_function *data = (base_function *) par;
	return _linear_interpolation(x, data->x, data->fx, data->points);
}

template<typename number>
number CustomInteraction<number>::_dfx(number x, void *par) {
	base_function *data = (base_function *) par;
	return _linear_interpolation(x, data->x, data->dfx, data->points);
}

template<typename number>
void CustomInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputString(&inp, "custom_lt_file", _lt_filename, 1);
	getInputInt(&inp, "custom_points", &_bonded_points, 0);
}

template<typename number>
void CustomInteraction<number>::init() {
	std::ifstream lt_file(_lt_filename, ios::in);
	if(!lt_file.good()) throw oxDNAException("Can't read lookup file '%s'. Aborting", _lt_filename);

	string line;
	int n_lines = 0;
	bool stop = false;
	while(!stop) {
		getline(lt_file, line);
		vector<string> spl = Utils::split(line);
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else n_lines++;
	}

	base_function data;
	data.x = new number[n_lines];
	data.fx = new number[n_lines];
	data.dfx = new number[n_lines];
	data.points = n_lines;

	lt_file.seekg(0, ios::beg);
	int i = 0;
	stop = false;
	while(!stop) {
		getline(lt_file, line);
		vector<string> spl = Utils::split(line);
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else {
			data.x[i] = atof(spl[0].c_str());
			data.fx[i] = atof(spl[1].c_str());
			data.dfx[i] = atof(spl[2].c_str());

			if(i > 0 && data.x[i] <= data.x[i-1]) throw oxDNAException("The x values of the lookup table should be monotonically increasing (found x[%d] = %f <= %f = x[%d])", i, data.x[i], i-1, data.x[i-1]);
			i++;
		}
	}
	lt_file.close();

	number lowlimit = data.x[0];
	number uplimit = data.x[i-1];

	this->_build_mesh(this, &CustomInteraction::_fx, &CustomInteraction::_dfx, (void *)(&data), _bonded_points, lowlimit, uplimit, _bonded_mesh);

//	ofstream out("yo.dat");
//	for(i = 0; i < 500; i++) {
//		number x = (uplimit - lowlimit)/500. * i;
//		number fx = this->_query_mesh(x, _bonded_mesh);
//		number dfx = this->_query_meshD(x, _bonded_mesh);
//		out << x << " " << fx << " " << dfx << endl;
//	}
//	out.close();

	delete[] data.x;
	delete[] data.fx;
	delete[] data.dfx;

	this->_rcut = uplimit;
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void CustomInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new CustomParticle<number>();
}

template<typename number>
void CustomInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) particles[i]->index = particles[i]->strand_id = i;

	CustomParticle<number> *p = (CustomParticle<number> *) particles[0];
	p->add_bonded_neigh((CustomParticle<number> *) particles[1]);
}

template<typename number>
number CustomInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number CustomInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = (number) 0.f;

	if(q == P_VIRTUAL) {
		CustomParticle<number> *Cp = (CustomParticle<number> *) p;
		for(typename set<CustomParticle<number> *>::iterator it = Cp->bonded_neighs.begin(); it != Cp->bonded_neighs.end(); it++) {
			energy += pair_interaction_bonded(p, *it, r, update_forces);
		}
	}
	else if(p->is_bonded(q)){
		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}

		if(p->index < q->index && update_forces) return energy;

		number dist = r->module();
		energy = this->_query_mesh(dist, _bonded_mesh);

		if(update_forces) {
			number force_mod = -this->_query_meshD(dist, _bonded_mesh);
			LR_vector<number> force = *r * (force_mod/dist);
			p->force -= force;
			q->force += force;
		}
	}

	return energy;
}

template<typename number>
number CustomInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return 0.f;
}

template<typename number>
void CustomInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class CustomInteraction<float>;
template class CustomInteraction<double>;

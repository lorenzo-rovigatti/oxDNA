/*
 * TSPAnalysis.cpp
 *
 *  Created on: 6 Nov 2014
 *      Author: lorenzo
 */

#include <sstream>

#include "TSPAnalysis.h"

using namespace std;

template<typename number>
TSPAnalysis<number>::TSPAnalysis(): BaseObservable<number>() {
	_N_stars = -1;
}

template<typename number>
TSPAnalysis<number>::~TSPAnalysis() {

}

template<typename number>
void TSPAnalysis<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputString(&sim_inp, "topology", _topology_filename, 1);
}

template<typename number>
void TSPAnalysis<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	ifstream topology(_topology_filename.c_str(), ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename.c_str());

	char line[512];
	topology.getline(line, 512);
	sscanf(line, "%d\n", &_N_stars);

	_stars.resize(_N_stars);
	for(int i = 0 ; i < *(config_info.N); i++) {
		TSPParticle<number> *p = (TSPParticle<number> *) config_info.particles[i];
		int ns = p->strand_id;

		if(p->is_anchor()) _stars[ns].set_anchor(p);
		else {
			if(p->arm() >= _stars[ns].n_arms()) _stars[ns].set_n_arms(p->arm()+1);
			_stars[ns].add_particle(p, p->arm());
		}
	}

	for(int i = 0; i < _N_stars; i++) _stars[i].set_config_info(this->_config_info);
}

template<typename number>
string TSPAnalysis<number>::get_output_string(llint curr_step) {
	stringstream ss;

	number avg_patches = 0.;

	for(int i = 0; i < _N_stars;i++) {
		_stars[i].update();
		avg_patches += _stars[i].patches.size();
	}
	avg_patches /= _N_stars;
	ss << avg_patches;

	return ss.str();
}

template class TSPAnalysis<float>;
template class TSPAnalysis<double>;

template<typename number>
void TSP<number>::set_n_arms(int na) {
	arms.resize(na);
}

template<typename number>
void TSP<number>::_flip_neighs(int arm, vector<vector<int> > &bond_map, vector<int> &arm_to_patch) {
	int na = n_arms();
	int p_a = arm_to_patch[arm];
	for(int i = 0; i < na; i++) {
		if(bond_map[arm][i] == 1) {
			int p_i = arm_to_patch[i];
			if(p_i > p_a) {
				arm_to_patch[i] = p_a;
				_flip_neighs(i, bond_map, arm_to_patch);
			}
		}
	}
}

template<typename number>
void TSP<number>::update() {
	int na = n_arms();

	vector<int> arm_to_patch(na);

	for(int i = 0; i < na; i++) {
		patches[i] = Patch<number>();
		arm_to_patch[i] = i;
	}

	vector<vector<int> > bond_map(na, vector<int>(na, 0));
	for(int i = 0; i < na; i++) {
		// c++ is retarded
		typename vector<TSPParticle<number> *>::iterator it;
		for(it = arms[i].begin(); it != arms[i].end(); it++) {
			TSPParticle<number> *p = *it;
			// we only care about attractive particles
			if(p->type != P_B) continue;
			vector<BaseParticle<number> *> neighs = _config_info.lists->get_neigh_list(p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				TSPParticle<number> *q = (TSPParticle<number> *) neighs[n];
				if(q->type == P_B && !p->is_bonded(q) && (p->arm() != q->arm())) {
					number energy = _config_info.interaction->pair_interaction_nonbonded(p, q, NULL, true);
					if(energy < 0.) bond_map[p->arm()][q->arm()] = bond_map[q->arm()][p->arm()] = 1;
				}
			}
		}
	}

	for(int i = 0; i < na; i++) _flip_neighs(i, bond_map, arm_to_patch);

	patches.clear();

	for(int i = 0; i < na; i++) {
		int patch = arm_to_patch[i];
		typename vector<TSPParticle<number> *>::iterator it;
		for(it = arms[i].begin(); it != arms[i].end(); it++) {
			TSPParticle<number> *p = *it;
			if(p->type == P_B) patches[patch].add_particle(p);
		}
	}

	for(int i = 0; i < na; i++) {
		patches[i].done();
		if(patches[i].empty()) patches.erase(i);
	}
}

template class TSP<float>;
template class TSP<double>;

template<typename number>
Patch<number>::Patch() : _pos(0., 0., 0.) {

}

template<typename number>
void Patch<number>::add_particle(TSPParticle<number> *p) {
	_particles.push_back(p);
	_arms.insert(p->arm());
}

template<typename number>
bool Patch<number>::empty() {
	return (_arms.size() < 2);
}

template<typename number>
void Patch<number>::done() {
	if(empty()) return;

	_arms.clear();
	_pos = LR_vector<number>(0., 0., 0.);

	typename vector<TSPParticle<number> *>::iterator it;
	for(it = _particles.begin(); it != _particles.end(); it++) {
		_pos += (*it)->pos;
	}

	_pos /= _particles.size();
}

template class Patch<float>;
template class Patch<double>;

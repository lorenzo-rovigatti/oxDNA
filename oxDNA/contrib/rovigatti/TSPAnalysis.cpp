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

	_attractive_monomers.resize(_N_stars);
	for(int i = 0 ; i < *config_info.N; i++) {
		BaseParticle<number> *p = config_info.particles[i];
		if(p->type == P_B) _attractive_monomers[i].push_back(p);
	}
}

template<typename number>
string TSPAnalysis<number>::get_output_string(llint curr_step) {
	stringstream ss;



	return ss.str();
}

template class TSPAnalysis<float>;
template class TSPAnalysis<double>;

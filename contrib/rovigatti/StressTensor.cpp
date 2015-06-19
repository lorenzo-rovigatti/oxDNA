/*
 * StressTensor.cpp
 *
 *  Created on: 11/ott/2013
 *      Author: lorenzo
 */

#include <sstream>

#include "StressTensor.h"
#include "FSInteraction.h"

using namespace std;

template<typename number>
StressTensor<number>::StressTensor() : BaseObservable<number>() {
	_print_averaged_off_diagonal = false;
}

template<typename number>
StressTensor<number>::~StressTensor() {

}

template<typename number>
void StressTensor<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	int tmp = 0;
	getInputBoolAsInt(&my_inp, "avg_off_diagonal", &tmp, 0);
	_print_averaged_off_diagonal = (bool)tmp;

	string topology_file;
	getInputString(&sim_inp, "topology", topology_file, 0);
	std::ifstream topology(topology_file.c_str(), ios::in);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%d %d\n", &_N, &_N_A);
	_N_B = _N - _N_A;
}

template<typename number>
string StressTensor<number>::get_output_string(llint curr_step) {
	vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();
	vector<vector<number> > st_pot(3, vector<number>(3, (number) 0));
	vector<vector<number> > st_kin(3, vector<number>(3, (number) 0));

	// we first set forces and torques to 0
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		p->force = LR_vector<number>(0, 0, 0);
		p->torque = LR_vector<number>(0, 0, 0);
	}

	// we reset and then update the 3-body lists, if present
	this->_config_info.interaction->pair_interaction_bonded(this->_config_info.particles[0], P_VIRTUAL);
	typename vector<ParticlePair<number> >::iterator it;
	for (it = pairs.begin(); it != pairs.end(); it ++ ) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;

		LR_vector<number> r = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
		this->_config_info.interaction->pair_interaction(p, q, &r, true);
	}

	// and then we compute the forces
	for (it = pairs.begin(); it != pairs.end(); it ++ ) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;

		LR_vector<number> old_p_force(p->force);
		LR_vector<number> old_q_force(q->force);
		LR_vector<number> old_p_torque(p->torque);
		LR_vector<number> old_q_torque(q->torque);

		p->force = q->force = p->torque = q->torque = LR_vector<number>();

		LR_vector<number> r = q->pos.minimum_image(p->pos, *this->_config_info.box_side);
		this->_config_info.interaction->pair_interaction(p, q, &r, true);
		number force_mod = q->force.module();
		number r_mod = r.module();

		// loop over all the 9 elements of the stress tensor
		for(int si = 0; si < 3; si++) {
			for(int sj = 0; sj < 3; sj++) {
				//st_pot[si][sj] += r[si] * q->force[sj];
				st_pot[si][sj] += r[si]*r[sj]*force_mod/r_mod;
			}
		}

		p->force = old_p_force;
		q->force = old_q_force;
		p->torque = old_p_torque;
		q->torque = old_q_torque;
	}

	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		for(int si = 0; si < 3; si++) {
			for(int sj = 0; sj < 3; sj++) {
				st_kin[si][sj] += p->vel[si]*p->vel[sj];
			}
		}
	}

	stringstream res;
	number mybox = *this->_config_info.box_side;
	int tot_N = 5*_N_A + 3*_N_B;

	res << curr_step << " " << curr_step << " " << tot_N << " " << 5*_N_A << " " << 0 << endl;
	res << mybox << " " << mybox << " " << mybox << " " << 0. << " " << 0. << " " << 0. << endl;

	//printf("%lf %lf\n", st_pot[0][0], other_st[0][0]);

	if(_print_averaged_off_diagonal) {
		number avg = 0;
		for (int si = 0; si < 3; si++) {
			for (int sj = 0; sj < 3; sj++) {
				if(si != sj) avg += st_pot[si][sj] + st_kin[si][sj];
			}
		}
		avg /= 6;
		res << avg;
	}
	else {
		FSInteraction<number> *fs_int = dynamic_cast<FSInteraction<number> *>(this->_config_info.interaction);
		if(fs_int != NULL) st_pot = fs_int->get_stress_tensor();

		number st_pot_01 = st_pot[0][1];
		number st_pot_02 = st_pot[0][2];
		number st_pot_12 = st_pot[1][2];
		/*number st_pot_01 = (st_pot[0][1] + st_pot[1][0])*0.5;
		number st_pot_02 = (st_pot[0][2] + st_pot[2][0])*0.5;
		number st_pot_12 = (st_pot[1][2] + st_pot[2][1])*0.5;*/
		res << st_pot[0][0] << " " << st_pot[1][1] << " " << st_pot[2][2] << " " << st_pot_01 << " " << st_pot_02 << " " << st_pot_12 << endl;

		number st_kin_01 = st_kin[0][1];
		number st_kin_02 = st_kin[0][2];
		number st_kin_12 = st_kin[1][2];
		/*number st_kin_01 = (st_kin[0][1] + st_kin[1][0])*0.5;
		number st_kin_02 = (st_kin[0][2] + st_kin[2][0])*0.5;
		number st_kin_12 = (st_kin[1][2] + st_kin[2][1])*0.5;*/
		res << st_kin[0][0] << " " << st_kin[1][1] << " " << st_kin[2][2] << " " << st_kin_01 << " " << st_kin_02 << " " << st_kin_12;
	}

	return res.str();
}

template class StressTensor<float>;
template class StressTensor<double>;

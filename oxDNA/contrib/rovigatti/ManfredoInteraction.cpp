/*
 * ManfredoInteraction.cpp
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#include "ManfredoInteraction.h"
#include <string>

using namespace std;

template<typename number>
ManfredoInteraction<number>::ManfredoInteraction() : BaseInteraction<number, ManfredoInteraction<number> >(), _N_per_tetramer(29) {
	//this->_int_map[LENNARD_JONES] = &ManfredoInteraction<number>::_lennard_jones;
	_N_tetramers = 0;
}

template<typename number>
ManfredoInteraction<number>::~ManfredoInteraction() {

}

template<typename number>
void ManfredoInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	_DNA_inter.get_settings(inp);

	getInputString(&inp, "lt_centre_arm_intra", _intra_filename[CENTRE_ARM_INTRA], 1);
	getInputInt(&inp, "lt_centre_arm_intra_points", &_intra_points[CENTRE_ARM_INTRA], 1);

	getInputString(&inp, "lt_arm_arm_near", _intra_filename[ARM_ARM_NEAR], 1);
	getInputInt(&inp, "lt_arm_arm_near_points", &_intra_points[ARM_ARM_NEAR], 1);

	getInputString(&inp, "lt_arm_arm_far", _intra_filename[ARM_ARM_FAR], 1);
	getInputInt(&inp, "lt_arm_arm_far_points", &_intra_points[ARM_ARM_FAR], 1);
}

template<typename number>
number ManfredoInteraction<number>::_linear_interpolation(number x, number *x_data, number *fx_data, int points) {
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
number ManfredoInteraction<number>::_fx(number x, void *par) {
	lt_data *data = (lt_data *) par;
	return _linear_interpolation(x, data->x, data->fx, data->points);
}

template<typename number>
number ManfredoInteraction<number>::_dfx(number x, void *par) {
	lt_data *data = (lt_data *) par;
	return _linear_interpolation(x, data->x, data->dfx, data->points);
}

template<typename number>
void ManfredoInteraction<number>::_build_intra_lt(int lt_type) {
	std::ifstream lt_file(_intra_filename[lt_type], ios::in);
	if(!lt_file.good()) throw oxDNAException("Can't read lookup file '%s'. Aborting", _intra_filename[lt_type]);

	string line;
	int n_lines = 0;
	bool stop = false;
	while(!stop) {
		getline(lt_file, line);
		vector<string> spl = Utils::split(line);
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else n_lines++;
	}

	lt_data data;
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

	this->_build_mesh(this, &ManfredoInteraction::_fx, &ManfredoInteraction::_dfx, (void *)(&data), _intra_points[lt_type], lowlimit, uplimit, _intra_mesh[lt_type]);

	delete[] data.x;
	delete[] data.fx;
	delete[] data.dfx;
}

template<typename number>
void ManfredoInteraction<number>::init() {
	_DNA_inter.init();

	_build_intra_lt(CENTRE_ARM_INTRA);
	_build_intra_lt(ARM_ARM_FAR);
	_build_intra_lt(ARM_ARM_NEAR);

	for(double x = 0; x < 25; x += 0.1) {
		printf("%lf %lf\n", x, _query_mesh(x, _intra_mesh[ARM_ARM_NEAR]));
	}
	exit(1);
}

template<typename number>
void ManfredoInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	int ind = 0;
	// for each tetramer
	for(int i = 0; i < _N_tetramers; i++) {
		particles[ind++] = new CustomParticle<number>();
		// for each arm
		for(int j = 0; j < 4; j++) {
			particles[ind++] = new CustomArmParticle<number>();
			//for(int k = 0; k < 6; k++) particles[ind++] = new DNANucleotide<number>(false);
		}
	}
}

template<typename number>
void ManfredoInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_tetramers);
	*N_strands = _N_tetramers;
	if(_N_tetramers*_N_per_tetramer != N) throw oxDNAException("Incoherent topology: the total number of particles should be equal to the number of tetramers (%d) times the number of particles in a tetramer (%d)", _N_tetramers, _N_per_tetramer);

	allocate_particles(particles, N);
	char arm_types[8] = "ACGATCG";
	for(int i = 0; i < N; i ++) {
		BaseParticle<number> *p = particles[i];
		p->index = i;
		p->strand_id = i / _N_per_tetramer;
		int p_type = _get_p_type(p);
		if(p_type == CENTRE) p->type = p->btype = P_A;
		else if(p_type == ARM) {
			// if it's in the arm then we assign it a letter, corresponding to a base type
			int rel_ind = p->index % _N_per_tetramer;
			int rel_arm = (rel_ind-1) % 7;
			p->type = p->btype = Utils::decode_base(arm_types[rel_arm]);
			// and then we assign 3" and 5" neighbours
			p->n3 = p-1;
			p->n5 = p+1;
			// if it's the first one it doesn't have the 3" neighbour
			if(rel_arm == 0) p->n3 = P_VIRTUAL;
			// and if it's the last one it doesn't have the 5"
			else if(rel_arm == 6) p->n5 = P_VIRTUAL;
		}
	}

	// here, we set the bonded neighbours inside each tetramer. The centre is bonded to each arm
	// and arms are bonded with each other
	for(int i = 0; i < _N_tetramers; i++) {
		CustomParticle<number> *centre = (CustomParticle<number> *) particles[i*_N_per_tetramer];
		for(int j = 0; j < 4; j++) {
			CustomParticle<number> *arm_j = (CustomParticle<number> *) particles[centre->index + 1 + 7*j];
			centre->add_bonded_neigh(arm_j);
			for(int k = j+1; k < 4; k++) {
				CustomParticle<number> *arm_k = (CustomParticle<number> *) particles[centre->index + 1 + 7*k];
				arm_j->add_bonded_neigh(arm_k);
			}
		}
	}
}

template<typename number>
int ManfredoInteraction<number>::_get_p_type(BaseParticle<number> *p) {
	int rel_ind = p->index % _N_per_tetramer;
	if(rel_ind == 0) return CENTRE;
	if(((rel_ind-1) % 7) == 0) return ARM;
	return STICKY;
}

template<typename number>
int ManfredoInteraction<number>::_get_inter_type(BaseParticle<number> *p, BaseParticle<number> *q) {
	int p_type = _get_p_type(p);
	int q_type = _get_p_type(p);

	// common for both inter- and intra-tetramer interactions
	if(_any(p_type, q_type, CENTRE) && _any(p_type, q_type, STICKY)) return CENTRE_STICKY;
	if(_any(p_type, q_type, ARM) && _any(p_type, q_type, STICKY)) return ARM_STICKY;
	if(_both(p_type, q_type, STICKY)) return STICKY_STICKY;

	// intra-tetramer interactions
	if(p->strand_id == q->strand_id) {
		if(_any(p_type, q_type, CENTRE) && _any(p_type, q_type, ARM)) return CENTRE_ARM_INTRA;
		if(_both(p_type, q_type, ARM)) {
			int diff = abs(p->index - q->index);
			if(diff == 14) return ARM_ARM_FAR;
			else return ARM_ARM_NEAR;
		}
	}
	// inter-tetramer interactions
	else {
		if(_any(p_type, q_type, CENTRE) && _any(p_type, q_type, ARM)) return CENTRE_ARM_INTER;
		if(_any(p_type, q_type, ARM) && _any(p_type, q_type, ARM)) return ARM_ARM_INTER;
		if(_both(p_type, q_type, CENTRE)) return CENTRE_CENTRE;
	}

	throw oxDNAException("No way we don't know which kind of interaction this is! We got p_type == %d and q_type == %d", p_type, q_type);
}

// if x is within the range [xlow, xupp] then this method calls BaseInteraction's method, otherwise
// it performs a linear extrapolation by using the nearest f(x) value (given by either A[0] or A[size-1]) and then
// by using the nearest f'(x) (given by B[0] or B[N-1]) as the slope for the linear interpolation.
template<typename number>
inline number ManfredoInteraction<number>::_query_mesh(number x, Mesh<number> &m) {
	if (x <= m.xlow) return m.A[0] + m.B[0]*(x - m.xlow);
	if (x >= m.xupp) return m.A[m.N-1] + m.B[m.N-1]*(x - m.xupp);
	return BaseInteraction<number, ManfredoInteraction<number> >::_query_mesh(x, m);
}

template<typename number>
number ManfredoInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number ManfredoInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	int p_type = _get_p_type(p);
	if(p_type == STICKY) return _DNA_inter.pair_interaction_bonded(p, q, r, update_forces);

	number energy = 0.f;
	if(q == P_VIRTUAL) {
		CustomParticle<number> *cp = (CustomParticle<number> *) p;
		for(typename set<CustomParticle<number> *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
			if(p->index > (*it)->index) energy += pair_interaction_bonded(p, *it, r, update_forces);
		}
	}
	else if(p->is_bonded(q)) {
		// bonded interactions are only intra-tetramer
		if(p->strand_id != q->strand_id) return energy;
		int type = _get_inter_type(p, q);

		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				// no periodic boundary conditions here, by choice
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}

		number dist = r->module();
		energy = this->_query_mesh(dist, _intra_mesh[type]);

		if(update_forces) {
			number force_mod = -this->_query_meshD(dist, _intra_mesh[type]);
			LR_vector<number> force = *r * (force_mod/dist);
			p->force -= force;
			q->force += force;
		}
	}
	return energy;
}

template<typename number>
number ManfredoInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	int p_type = _get_p_type(p);
	int q_type = _get_p_type(q);
	if(_both(p_type, q_type, STICKY)) return _DNA_inter.pair_interaction_nonbonded(p, q, r, update_forces);

	return 0.f;
}

template<typename number>
void ManfredoInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

// we copy the next method from the DNANucleotide class
// not a pretty way of doing it, but it works. Sorry.
template<typename number>
void CustomArmParticle<number>::set_positions() {
	this->int_centers[DNANucleotide<number>::BACK] = this->orientation*_principal_DNA_axis*POS_BACK;
	this->int_centers[DNANucleotide<number>::STACK] = this->int_centers[DNANucleotide<number>::BACK]*(POS_STACK/POS_BACK);
	this->int_centers[DNANucleotide<number>::BASE] = this->int_centers[DNANucleotide<number>::STACK]*(POS_BASE/POS_STACK);
}

template class CustomArmParticle<float>;
template class CustomArmParticle<double>;

template class ManfredoInteraction<float>;
template class ManfredoInteraction<double>;

/*
 * ManfredoInteraction.cpp
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#include "ManfredoInteraction.h"
#include <string>

using namespace std;

ManfredoInteraction::ManfredoInteraction() :
				BaseInteraction<ManfredoInteraction>(),
				_N_per_tetramer(29) {
	_N_tetramers = 0;
	_int_map[0] = &ManfredoInteraction::pair_interaction;
}

ManfredoInteraction::~ManfredoInteraction() {

}

void ManfredoInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	_DNA_inter.get_settings(inp);

	char my_T[512];
	getInputString(&inp, "T", my_T, 1);
	_T = Utils::get_temperature(my_T);

	getInputString(&inp, "lt_centre_arm_intra", _intra_filename[CENTRE_ARM_INTRA], 1);
	getInputInt(&inp, "lt_centre_arm_intra_points", &_intra_points[CENTRE_ARM_INTRA], 1);

	getInputString(&inp, "lt_arm_arm_near", _intra_filename[ARM_ARM_NEAR], 1);
	getInputInt(&inp, "lt_arm_arm_near_points", &_intra_points[ARM_ARM_NEAR], 1);

	getInputString(&inp, "lt_arm_arm_far", _intra_filename[ARM_ARM_FAR], 1);
	getInputInt(&inp, "lt_arm_arm_far_points", &_intra_points[ARM_ARM_FAR], 1);

	getInputString(&inp, "lt_centre_centre", _inter_filename[CENTRE_CENTRE], 1);
	getInputInt(&inp, "lt_centre_centre_points", &_inter_points[CENTRE_CENTRE], 1);

	getInputString(&inp, "lt_centre_arm", _inter_filename[CENTRE_ARM_INTER], 1);
	getInputInt(&inp, "lt_centre_arm_points", &_inter_points[CENTRE_ARM_INTER], 1);

	getInputString(&inp, "lt_arm_arm", _inter_filename[ARM_ARM_INTER], 1);
	getInputInt(&inp, "lt_arm_arm_points", &_inter_points[ARM_ARM_INTER], 1);

	_rcut = 25.;
}

number ManfredoInteraction::_linear_interpolation(number x, number *x_data, number *fx_data, int points) {
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

number ManfredoInteraction::_fx(number x, void *par) {
	lt_data *data = (lt_data *) par;
	return _linear_interpolation(x, data->x, data->fx, data->points);
}

number ManfredoInteraction::_dfx(number x, void *par) {
	lt_data *data = (lt_data *) par;
	return _linear_interpolation(x, data->x, data->dfx, data->points);
}

void ManfredoInteraction::_build_lt(Mesh &mesh, int points, char *filename) {
	std::ifstream lt_file(filename, ios::in);
	if(!lt_file.good()) throw oxDNAException("Can't read lookup file '%s'. Aborting", filename);

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

	lt_file.clear();
	lt_file.seekg(0, ios::beg);
	int i = 0;
	stop = false;
	while(!stop) {
		getline(lt_file, line);
		vector<string> spl = Utils::split(line);
		if(spl.size() != 3 || lt_file.eof()) stop = true;
		else {
			data.x[i] = atof(spl[0].c_str());
			data.fx[i] = _T * atof(spl[1].c_str());
			data.dfx[i] = -_T * atof(spl[2].c_str());

			if(i > 0 && data.x[i] <= data.x[i - 1]) throw oxDNAException("The x values of the lookup table should be monotonically increasing (found x[%d] = %f <= %f = x[%d])", i, data.x[i], i - 1, data.x[i - 1]);
			i++;
		}
	}
	lt_file.close();

	number lowlimit = data.x[0];
	number uplimit = data.x[i - 1];

	_build_mesh(this, &ManfredoInteraction::_fx, &ManfredoInteraction::_dfx, (void *) (&data), points, lowlimit, uplimit, mesh);

	delete[] data.x;
	delete[] data.fx;
	delete[] data.dfx;
}

void ManfredoInteraction::init() {
	_DNA_inter.init();

	for(int i = 0; i < 3; i++)
		_build_lt(_intra_mesh[i], _intra_points[i], _intra_filename[i]);
	for(int i = 0; i < 3; i++)
		_build_lt(_inter_mesh[i], _inter_points[i], _inter_filename[i]);
}

void ManfredoInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	int ind = 0;
	// for each tetramer
	for(int i = 0; i < _N_tetramers; i++) {
		particles[ind++] = new CustomParticle();
		// for each arm
		for(int j = 0; j < 4; j++) {
			particles[ind++] = new CustomArmParticle();
			for(int k = 0; k < 6; k++)
				particles[ind++] = new DNANucleotide(false);
		}
	}
}

void ManfredoInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	sscanf(line, "%*d %d\n", &_N_tetramers);
	*N_strands = _N_tetramers;
	if(_N_tetramers * _N_per_tetramer != N) {
		throw oxDNAException("Incoherent topology: the total number of particles should be equal to the number of tetramers (%d) times the number of particles in a tetramer (%d)", _N_tetramers, _N_per_tetramer);
	}

	allocate_particles(particles);
	char arm_types[8] = "ACGATCG";
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		p->index = i;
		p->strand_id = i / _N_per_tetramer;
		int p_type = _get_p_type(p);
		if(p_type == CENTRE) p->type = p->btype = P_A;
		else {
			// if it's in the arm then we assign it a letter, corresponding to a base type
			int rel_ind = p->index % _N_per_tetramer;
			int rel_arm = (rel_ind - 1) % 7;
			p->type = p->btype = Utils::decode_base(arm_types[rel_arm]);
			// and then we assign 3" and 5" neighbours
			p->n3 = particles[i - 1];
			p->n5 = particles[i + 1];
			// if it's the first one it doesn't have the 3" neighbour
			if(rel_arm == 0) p->n3 = P_VIRTUAL;
			// and if it's the last one it doesn't have the 5"
			else if(rel_arm == 6) p->n5 = P_VIRTUAL;
		}
	}

	// here, we set the bonded neighbours inside each tetramer. The centre is bonded to each arm
	// and arms are bonded with each other
	for(int i = 0; i < _N_tetramers; i++) {
		CustomParticle *centre = (CustomParticle *) particles[i * _N_per_tetramer];
		for(int j = 0; j < 4; j++) {
			CustomParticle *arm_j = (CustomParticle *) particles[centre->index + 1 + 7 * j];
			centre->add_bonded_neigh(arm_j);
			for(int k = j + 1; k < 4; k++) {
				CustomParticle *arm_k = (CustomParticle *) particles[centre->index + 1 + 7 * k];
				arm_j->add_bonded_neigh(arm_k);
			}
		}
	}
}

int ManfredoInteraction::_get_p_type(BaseParticle *p) {
	int rel_ind = p->index % _N_per_tetramer;
	if(rel_ind == 0) return CENTRE;
	if(((rel_ind - 1) % 7) == 0) return ARM;

	return STICKY;
}

int ManfredoInteraction::_get_inter_type(BaseParticle *p, BaseParticle *q) {
	int p_type = _get_p_type(p);
	int q_type = _get_p_type(q);

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
		if(_both(p_type, q_type, ARM)) return ARM_ARM_INTER;
		if(_both(p_type, q_type, CENTRE)) return CENTRE_CENTRE;
	}

	throw oxDNAException("No way we don't know which kind of interaction this is! We got p_type == %d and q_type == %d", p_type, q_type);
}

// if x is within the range [xlow, xupp] then this method calls BaseInteraction's method, otherwise
// it performs a linear extrapolation by using the nearest f(x) value (given by either A[0] or A[size-1]) and then
// by using the nearest f'(x) (given by B[0] or B[N-1]) as the slope for the linear interpolation.

inline number ManfredoInteraction::_query_mesh(number x, Mesh &m) {
	if(x <= m._xlow) return m._A[0] + m._B[0] * (x - m._xlow);
	if(x >= m._xupp) return m._A[m.N - 1] + m._B[m.N - 1] * (x - m._xupp);
	return BaseInteraction<ManfredoInteraction>::_query_mesh(x, m);
}

inline number ManfredoInteraction::_query_meshD(number x, Mesh &m) {
	if(x <= m._xlow) return m._B[0];
	if(x >= m._xupp) return m._B[m.N - 1];
	return BaseInteraction<ManfredoInteraction>::_query_meshD(x, m);
}

number ManfredoInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number ManfredoInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.f;
	if(q == P_VIRTUAL) {
		// if the particle is a sticky end nucleotide then we have to invoke the DNA_inter method to handle all
		// its bonded interactions
		if(_get_p_type(p) == STICKY) return _DNA_inter.pair_interaction_bonded(p, q, compute_r, update_forces);

		CustomParticle *cp = (CustomParticle *) p;
		for(typename set<CustomParticle *>::iterator it = cp->bonded_neighs.begin(); it != cp->bonded_neighs.end(); it++) {
			if(p->index > (*it)->index) {
				energy += pair_interaction_bonded(p, *it, true, update_forces);
			}
		}
	}
	else if(p->is_bonded(q)) {
		// bonded interactions are only intra-tetramer
		if(p->strand_id != q->strand_id) return energy;

		int type = _get_inter_type(p, q);

		if(type == STICKY_STICKY || type == ARM_STICKY) {
			if(update_forces) {
				throw oxDNAException("STICKY_STICKY and ARM_STICKY should always be handled in the q == P_VIRTUAL portion\n");
			}
			return _DNA_inter.pair_interaction_bonded(p, q, compute_r, update_forces);
		}
		if(compute_r) {
			if(q != P_VIRTUAL && p != P_VIRTUAL) {
				// no periodic boundary conditions here, by choice
				_computed_r = q->pos - p->pos;
			}
		}

		number dist = _computed_r.module();
		energy = _query_mesh(dist, _intra_mesh[type]);

		if(update_forces) {
			number force_mod = -_query_meshD(dist, _intra_mesh[type]);
			LR_vector force = _computed_r * (force_mod / dist);
			p->force -= force;
			q->force += force;
		}
	}
	return energy;
}

number ManfredoInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.f;

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	int p_type = _get_p_type(p);
	int q_type = _get_p_type(q);

	int type = _get_inter_type(p, q);

	if(type == STICKY_STICKY || type == ARM_STICKY) {
		return _DNA_inter.pair_interaction_nonbonded(p, q, false, update_forces);
	}

	// if there are no sticky ends involved and if the two particles belong to two different tetramers
	if(!_any(p_type, q_type, STICKY) && p->strand_id != q->strand_id) {
		number dist = _computed_r.module();
		energy = _query_mesh(dist, _inter_mesh[type]);

		if(update_forces) {
			number force_mod = -_query_meshD(dist, _inter_mesh[type]);
			LR_vector force = _computed_r * (force_mod / dist);
			p->force -= force;
			q->force += force;
		}
	}

	return energy;
}

void ManfredoInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

void ManfredoInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	number sqr_limit = SQR(20.);
	LR_vector arm_poss[4] = { LR_vector(-3.66451, -2.22875, -7.98832), LR_vector(3.83406, 7.82627, 1.45587), LR_vector(-3.6474, 2.70588, 7.63758), LR_vector(3.41029, -7.53441, -3.37317) };

	LR_vector box_sides = _box->box_sides();
	std::vector<LR_vector> tetra_centres;
	int p_ind = 0;
	for(int nt = 0; nt < _N_tetramers; nt++) {
		BaseParticle *centre = particles[p_ind];
		bool found = false;
		while(!found) {
			found = true;
			centre->pos = Utils::get_random_vector();
			centre->pos.x *= box_sides.x;
			centre->pos.y *= box_sides.y;
			centre->pos.z *= box_sides.z;
			typename std::vector<LR_vector>::iterator it;
			for(it = tetra_centres.begin(); it != tetra_centres.end() && found; it++) {
				number sqr_dist = _box->sqr_min_image_distance(*it, centre->pos);
				if(sqr_dist < sqr_limit) found = false;
			}
		}
		tetra_centres.push_back(centre->pos);
		p_ind++;

		for(int na = 0; na < 4; na++) {
			BaseParticle *arm = particles[p_ind];
			arm->pos = centre->pos + arm_poss[na];
			//arm->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
			//arm->set_positions();

			p_ind++;
			LR_vector dir = arm_poss[na];
			dir.normalize();

			LR_matrix orientation(dir, Utils::get_random_vector(), Utils::get_random_vector());
			Utils::orthonormalize_matrix(orientation);
			if(orientation.determinant() < 0.) orientation.v1 = -orientation.v1;
			LR_vector tmp = orientation.v1;
			orientation.v1 = orientation.v3;
			orientation.v3 = orientation.v2;
			orientation.v2 = tmp;
			orientation.transpone();

			arm->orientation = orientation;
			arm->set_positions();

			for(int nsticky = 0; nsticky < 6; nsticky++) {
				BaseParticle *sticky = particles[p_ind];
				//sticky->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
				sticky->orientation = orientation;
				sticky->set_positions();

				LR_vector prev_back = particles[p_ind - 1]->pos + particles[p_ind - 1]->int_centers[DNANucleotide::BACK];
				LR_vector curr_back = prev_back + dir * FENE_R0_OXDNA;
				sticky->pos = curr_back - sticky->int_centers[DNANucleotide::BACK];
				p_ind++;
			}
		}
	}
}

CustomArmParticle::CustomArmParticle() :
				CustomParticle(),
				_principal_DNA_axis(LR_vector(1, 0, 0)) {
	int_centers.resize(3);
}
;

// we copy the next method from the DNANucleotide class
// not a pretty way of doing it, but it works. Sorry.

void CustomArmParticle::set_positions() {
	int_centers[DNANucleotide::BACK] = orientation * _principal_DNA_axis * POS_BACK;
	int_centers[DNANucleotide::STACK] = int_centers[DNANucleotide::BACK] * (POS_STACK / POS_BACK);
	int_centers[DNANucleotide::BASE] = int_centers[DNANucleotide::STACK] * (POS_BASE / POS_STACK);
}

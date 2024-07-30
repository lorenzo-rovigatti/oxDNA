/*
 * TSPAnalysis.cpp
 *
 *  Created on: 6 Nov 2014
 *      Author: lorenzo
 */

#include <sstream>

#include "TSPAnalysis.h"

using namespace std;

TSPAnalysis::TSPAnalysis() :
				BaseObservable() {
	_N_stars = -1;
	_mode = TSP_ALL;
	_is_SPB = false;
	_sqr_rg = _delta = _S = _b = _c = 0.;
	_t = 0;
}

TSPAnalysis::~TSPAnalysis() {

}

void TSPAnalysis::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	string int_type;
	getInputString(&sim_inp, "interaction_type", int_type, 0);
	if(int_type.compare("GraftedInteraction") == 0) _is_SPB = true;
	if(!_is_SPB) getInputString(&sim_inp, "topology", _topology_filename, 1);

	std::string my_mode;
	if(getInputString(&my_inp, "mode", my_mode, 0) == KEY_FOUND) {
		if(my_mode == "size") _mode = TSP_SP;
		else if(my_mode == "all") _mode = TSP_ALL;
		else if(my_mode == "all_angles") _mode = TSP_ALL_ANGLES;
		else if(my_mode == "eigenvalues") _mode = TSP_EIGENVALUES;
		else throw oxDNAException("TSPAnalysis: Mode '%s' not supported", my_mode.c_str());
	}
}

void TSPAnalysis::init() {
	BaseObservable::init();

	if(!_is_SPB) {
		ifstream topology(_topology_filename.c_str(), ios::in);
		if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename.c_str());

		char line[512];
		topology.getline(line, 512);
		sscanf(line, "%d\n", &_N_stars);
		topology.close();
	}
	else _N_stars = 1;

	_stars.resize(_N_stars, TSP(_is_SPB));
	for(int i = 0; i < _config_info->N(); i++) {
		TSPParticle *p = (TSPParticle *) _config_info->particles()[i];
		int ns = (!_is_SPB) ? p->strand_id : 0;
		if(ns >= _N_stars) throw oxDNAException("ns (%d) >= _N_stars (%d)", ns, _N_stars);

		if(p->is_anchor() || (_is_SPB && i == 0)) _stars[ns].set_anchor(p);
		else {
			if(p->arm() >= _stars[ns].n_arms()) _stars[ns].set_n_arms(p->arm() + 1);
			_stars[ns].add_particle(p, p->arm());
		}
	}

	for(int i = 0; i < _N_stars; i++) {
		_stars[i].set_config_info(_config_info);
	}
}

string TSPAnalysis::get_output_string(llint curr_step) {
	stringstream ss;
	_t++;

	number avg_patches = 0.;
	number avg_distance = 0.;
	number avg_patch_size = 0.;
	number avg_angle = 0.;

	for(int i = 0; i < _N_stars; i++) {
		_stars[i].update();
		if(_mode == TSP_EIGENVALUES) {
			number l1 = _stars[i].l1;
			number l2 = _stars[i].l2;
			number l3 = _stars[i].l3;

			ss << l1 << " " << l2 << " " << l3 << " ";

			number I1 = l1 + l2 + l3;
			number I2 = l1 * l2 + l2 * l3 + l3 * l1;
			//number I3 = l1*l2*l3;

			_sqr_rg += I1;
			_delta += 1. - 3 * I2 / SQR(I1);
			_S += (3 * l1 - I1) * (3 * l2 - I1) * (3 * l3 - I1) / CUB(I1);
			_b += l1 - 0.5 * (l2 + l3);
			_c += l2 - l3;

			number f = _N_stars * _t;
			number rg2 = _sqr_rg / f;

			ss << _delta / f << " " << _S / f << " " << _b / f / rg2 << " " << _c / f / rg2 << " " << sqrt(rg2) << endl;

			continue;
		}
		int n_patches = _stars[i].patches.size();
		avg_patches += n_patches;
		LR_vector star_pos = _stars[i].pos();
		number avg_star_distance = 0.;
		number avg_star_patch_size = 0.;
		number avg_star_angle = 0.;

		if(n_patches > 0) {
			typename map<int, Patch>::iterator it;
			for(it = _stars[i].patches.begin(); it != _stars[i].patches.end(); it++) {
				Patch &patch = it->second;
				LR_vector rel_pos_1 = _config_info->box->min_image(patch.pos(), star_pos);
				number dist_1 = sqrt(SQR(rel_pos_1));
				avg_star_distance += dist_1;
				avg_star_patch_size += patch.n_arms();

				if(_mode == TSP_SP) ss << patch.n_arms() << " " << dist_1 << endl;

				rel_pos_1 /= dist_1;
				// angle
				typename map<int, Patch>::iterator jt;
				for(jt = it; jt != _stars[i].patches.end(); jt++) {
					if(it != jt) {
						Patch &patch_2 = jt->second;
						LR_vector rel_pos_2 = _config_info->box->min_image(patch_2.pos(), star_pos);
						rel_pos_2 /= sqrt(SQR(rel_pos_2));
						number angle = acos(rel_pos_1 * rel_pos_2);

						if(_mode == TSP_ALL_ANGLES) ss << angle << " " << patch.n_arms() << " " << patch_2.n_arms() << endl;

						avg_star_angle += angle;
					}
				}
			}
			avg_distance += avg_star_distance / (number) n_patches;
			avg_patch_size += avg_star_patch_size / (number) n_patches;
			// the factor 2 is to count for the number of unique pairs, which is N*(N-1) / 2
			if(n_patches > 1) avg_angle += 2. * avg_star_angle / (number)(n_patches * (n_patches - 1));
		}
	}

	avg_patches /= _N_stars;
	avg_distance /= _N_stars;
	avg_patch_size /= _N_stars;
	avg_angle /= _N_stars;

	if(_mode == TSP_ALL) ss << avg_patches << " " << avg_patch_size << " " << avg_distance << " " << avg_angle;

	return ss.str();
}

void TSP::set_n_arms(int na) {
	arms.resize(na);
}

void TSP::_flip_neighs(int arm, vector<vector<int> > &bond_map, vector<int> &arm_to_patch) {
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

void TSP::_compute_eigenvalues() {
	LR_vector com(0, 0, 0);

	int na = n_arms();
	int N = 0;
	for(int i = 0; i < na; i++) {
		typename vector<TSPParticle *>::iterator it;
		for(it = arms[i].begin(); it != arms[i].end(); it++) {
			TSPParticle *p = *it;
			com += _config_info->box->get_abs_pos(p);
			N++;
		}
	}
	com /= N;

	double IM[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	for(int i = 0; i < na; i++) {
		typename vector<TSPParticle *>::iterator it;
		for(it = arms[i].begin(); it != arms[i].end(); it++) {
			TSPParticle *p = *it;
			LR_vector i_pos = _config_info->box->get_abs_pos(p) - com;

			IM[0][0] += SQR(i_pos[1]) + SQR(i_pos[2]);
			IM[0][1] += -i_pos[0] * i_pos[1];
			IM[0][2] += -i_pos[0] * i_pos[2];

			IM[1][1] += SQR(i_pos[0]) + SQR(i_pos[2]);
			IM[1][2] += -i_pos[1] * i_pos[2];

			IM[2][2] += SQR(i_pos[0]) + SQR(i_pos[1]);
		}
	}
	IM[1][0] = IM[0][1];
	IM[2][0] = IM[0][2];
	IM[2][1] = IM[1][2];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			IM[i][j] /= N;

	double EV[3][3];
	double val[3];
	eigen_decomposition(IM, EV, val);

	l1 = val[2];
	l2 = val[1];
	l3 = val[0];
}

void TSP::update() {
	_compute_eigenvalues();

	int na = n_arms();

	vector<int> arm_to_patch(na);

	for(int i = 0; i < na; i++) {
		patches[i] = Patch();
		arm_to_patch[i] = i;
	}

	vector<vector<int> > bond_map(na, vector<int>(na, 0));
	for(int i = 0; i < na; i++) {
		// c++ is retarded
		typename vector<TSPParticle *>::iterator it;
		for(it = arms[i].begin(); it != arms[i].end(); it++) {
			TSPParticle *p = *it;
			// we only care about attractive particles
			if(p->type != P_B) continue;
			vector<BaseParticle *> neighs = _config_info->lists->get_neigh_list(p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				TSPParticle *q = (TSPParticle *) neighs[n];
				// particles should be attractive, non bonded, on the same star but not on the same arm
				bool same_star = (!_is_SPB) ? (p->strand_id == q->strand_id) : true;
				if(q->type == P_B && !p->is_bonded(q) && same_star && (p->arm() != q->arm())) {
					number energy = _config_info->interaction->pair_interaction_nonbonded(p, q, true, true);
					if(energy < 0.) bond_map[p->arm()][q->arm()] = bond_map[q->arm()][p->arm()] = 1;
				}
			}
		}
	}

	for(int i = 0; i < na; i++)
		_flip_neighs(i, bond_map, arm_to_patch);

	patches.clear();

	for(int i = 0; i < na; i++) {
		int patch = arm_to_patch[i];
		typename vector<TSPParticle *>::iterator it;
		for(it = arms[i].begin(); it != arms[i].end(); it++) {
			TSPParticle *p = *it;
			if(p->type == P_B) patches[patch].add_particle(p);
		}
	}

	for(int i = 0; i < na; i++) {
		patches[i].done(_config_info->box);
		if(patches[i].empty()) patches.erase(i);
	}
}

Patch::Patch() :
				_pos(0., 0., 0.) {

}

void Patch::add_particle(TSPParticle *p) {
	_particles.push_back(p);
	_arms.insert(p->arm());
}

bool Patch::empty() {
	return (_arms.size() < 2);
}

void Patch::done(BaseBox *box) {
	if(empty()) return;

	_pos = LR_vector(0., 0., 0.);

	typename vector<TSPParticle *>::iterator it;
	for(it = _particles.begin(); it != _particles.end(); it++) {
		_pos += box->get_abs_pos(*it);
	}

	_pos /= _particles.size();
}

// STUFF TO DIAGONALIZE 3X3 MATRICES
#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static double hypot2(double x, double y) {
	return sqrt(x * x + y * y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

	for(int j = 0; j < n; j++) {
		d[j] = V[n - 1][j];
	}

	// Householder reduction to tridiagonal form.

	for(int i = n - 1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for(int k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if(scale == 0.0) {
			e[i] = d[i - 1];
			for(int j = 0; j < i; j++) {
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		}
		else {

			// Generate Householder vector.

			for(int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i - 1];
			double g = sqrt(h);
			if(f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i - 1] = f - g;
			for(int j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for(int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for(int k = j + 1; k <= i - 1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for(int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for(int j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for(int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for(int k = j; k <= i - 1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.

	for(int i = 0; i < n - 1; i++) {
		V[n - 1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i + 1];
		if(h != 0.0) {
			for(int k = 0; k <= i; k++) {
				d[k] = V[k][i + 1] / h;
			}
			for(int j = 0; j <= i; j++) {
				double g = 0.0;
				for(int k = 0; k <= i; k++) {
					g += V[k][i + 1] * V[k][j];
				}
				for(int k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for(int k = 0; k <= i; k++) {
			V[k][i + 1] = 0.0;
		}
	}
	for(int j = 0; j < n; j++) {
		d[j] = V[n - 1][j];
		V[n - 1][j] = 0.0;
	}
	V[n - 1][n - 1] = 1.0;
	e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

	for(int i = 1; i < n; i++) {
		e[i - 1] = e[i];
	}
	e[n - 1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	double eps = pow(2.0, -52.0);
	for(int l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = MAX(tst1, fabs(d[l]) + fabs(e[l]));
		int m = l;
		while(m < n) {
			if(fabs(e[m]) <= eps * tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if(m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l + 1] - g) / (2.0 * e[l]);
				double r = hypot2(p, 1.0);
				if(p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l + 1] = e[l] * (p + r);
				double dl1 = d[l + 1];
				double h = g - d[l];
				for(int i = l + 2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l + 1];
				double s = 0.0;
				double s2 = 0.0;
				for(int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p, e[i]);
					e[i + 1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i + 1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for(int k = 0; k < n; k++) {
						h = V[k][i + 1];
						V[k][i + 1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while(fabs(e[l]) > eps * tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.

	for(int i = 0; i < n - 1; i++) {
		int k = i;
		double p = d[i];
		for(int j = i + 1; j < n; j++) {
			if(d[j] < p) {
				k = j;
				p = d[j];
			}
		}
		if(k != i) {
			d[k] = d[i];
			d[i] = p;
			for(int j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}

void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
	double e[n];
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			V[i][j] = A[i][j];
		}
	}
	tred2(V, d, e);
	tql2(V, d, e);
}

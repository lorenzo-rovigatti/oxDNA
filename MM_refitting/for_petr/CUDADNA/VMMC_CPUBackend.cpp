/*
 *  Created on: 26/nov/2010
 *      Author: flavio
 */

#include "VMMC_CPUBackend.h"
#include "IOManager.h"
#include "Utils.h"
#include <set>

inline int find_i(int * clust, int size, int value) {
	int i;
	for (i =0; i < size; i ++) {
		if (clust[i] == value) {
			return i;
		}
	}
	return -1;
}

template<typename number> VMMC_CPUBackend<number>::VMMC_CPUBackend(IOManager *IO) : MC_CPUBackend<number> (IO) {
	//_op = NULL;
	_have_us = false;
	_netemps = 0;
	_etemps = NULL;
	_maxclust = 0;
	_reject_prelinks = false;
	_preserve_topology = false;
	_small_system = false;
	_last_move = MC_MOVE_TRANSLATION;
}

template<typename number>
VMMC_CPUBackend<number>::~VMMC_CPUBackend() {
	
	// this is because otherwise the pointer to the force object gets freed
	// twice... maybe not the best way to do this, but oh well...
	for (int i = 0; i < this->_N; i++) {
		_particles_old[i]._N_ext_forces = 0;
	}

	if (_particles_old != NULL) 
		delete[] _particles_old;
	
	_delete_cells();

	delete[] new_en3s;
	delete[] new_en5s;
	delete[] new_stn3s;
	delete[] new_stn5s;

	if (this->_N_updates > 0)
		divide_given_timing(&this->_timer, 1, this->_N / (double) this->_N_updates);
	
	if (_netemps > 0)
		delete[] _etemps;
	
	if (_small_system) {
		for (int k = 0; k < this->_N; k ++) {
			delete[] eijm[k];
			delete[] eijm_old[k];
			delete[] hbijm[k];
			delete[] hbijm_old[k];
		}
		delete[] eijm;
		delete[] eijm_old;
		delete[] hbijm;
		delete[] hbijm_old;
	}
	return;
}

template<typename number>
//void VMMC_CPUBackend<number>::init(ifstream &conf_input) {
void VMMC_CPUBackend<number>::init(char conf_filename[256]) {
	//MCBackend<number>::init(conf_input);
	MCBackend<number>::init (conf_filename);

	// fix maxclust if evidently wrong
	if (_maxclust < 1) {
		this->_IO->log(this->_IO->LOG_WARNING, "maxclust < 0, setting it to N = %i", this->_N);
		_maxclust = this->_N; 
	}
	if (_maxclust > this->_N) {
		this->_IO->log(this->_IO->LOG_WARNING, "maxclust > N does not make sense, setting it to N = %i", this->_N);
		_maxclust = this->_N;
	}

	if (_have_us) {
		_op.init_from_file(_op_file, this->_particles, this->_N, this->_IO);
		_w.init((const char *) _weights_file, &_op);
		if (_reload_hist)
			_h.init(_init_hist_file, &_op, _etemps, _netemps);
		else
			_h.init(&_op, _etemps, _netemps);
		_h.set_simtemp(this->_T);
	}

	for (int i = 0; i < this->_N; i++) {
		// this is needed for the first _compute_energy()
		this->_particles[i].set_positions();
		this->_particles[i].orientationT
				= this->_particles[i].orientation.get_transpose();
	}

	if (this->_delta[MC_MOVE_TRANSLATION] * sqrt(3) > this->_verlet_skin)
		this->_IO->die(
				"verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");
	
	// setting the maximum displacement
	if (_preserve_topology) {
		_max_move_size = (number) 0.5;
		_max_move_size_sqr = _max_move_size * _max_move_size;
		this->_IO->log(this->_IO->LOG_INFO, "Preserving topology; max_move_size = %lf...", _max_move_size);
	}
	else {
		_max_move_size = this->_box_side / 2. - 2.;
		_max_move_size_sqr = _max_move_size * _max_move_size;
		this->_IO->log(this->_IO->LOG_INFO, "Not attempting to preserve topology; max_move_size = %g", _max_move_size);
	}
	
	_particles_old = new Particle<number>[this->_N];
	for(int i = 0; i < this->_N; i++) {
		_particles_old[i].index = i;
		_particles_old[i].type = 1;
		_particles_old[i].init (this->_max_neigh);
		for (int l = 0; l < this->_particles[i]._N_ext_forces; l ++) {
			_particles_old[i].add_ext_force (this->_particles[i]._ext_forces[l]);
		}
		_particles_old[i].copy_from (this->_particles[i]);
		// set the initial forces
	}

	new_en3s = new number[this->_N];
	new_en5s = new number[this->_N];
	new_stn3s = new number[this->_N];
	new_stn5s = new number[this->_N];
	for (int k = 0; k < this->_N; k ++) {
		new_en3s[k] = new_en5s[k] = (number)0.; 
		new_stn3s[k] = new_stn5s[k] = (number)0.;
	}
	number tmpf, epq;
	Particle<number> * p, *q;
	for (int k = 0; k < this->_N; k ++) {
		p = &this->_particles[k];
		if (p->n3 != P_VIRTUAL) {
			q = &this->_particles[p->n3];
			epq = _particle_particle_bonded_interaction_n3 (p, q, &tmpf);
			p->en3 = epq;
			q->en5 = epq;
			p->esn3 = tmpf;
			q->esn5 = tmpf;
		}
	}

	
	if (_small_system) {
		eijm = new number*[this->_N];
		eijm_old = new number*[this->_N];
		hbijm = new bool*[this->_N];
		hbijm_old = new bool*[this->_N];
		for (int k = 0; k < this->_N; k ++) {
			eijm[k] = new number[this->_N];
			eijm_old[k] = new number[this->_N];
			hbijm[k] = new bool[this->_N];
			hbijm_old[k] = new bool[this->_N];
		}
		/*eijm = new number[this->_N][this->_N];
		eijm_old = new number[this->_N][this->_N];
		hbijm = new bool[this->_N][this->_N];
		hbijm_old = new bool[this->_N][this->_N];
		*/
		
		for (int k = 0; k < this->_N; k ++) {
			for (int l = 0; l < k; l ++) {
				p = &this->_particles[k];
				q = &this->_particles[l];
				if (p->n3 != q->index && p->n5 != q->index) {
					eijm[k][l] = eijm[l][k] = eijm_old[k][l] = eijm_old[l][k] = _particle_particle_interaction (p, q, &tmpf);
					hbijm[k][l] = hbijm[l][k] = hbijm_old[k][l] = hbijm_old[l][k] = (tmpf < HB_CUTOFF);
				}
			}
		}
	}

	_init_cells();
	_create_cells();
	
	this->_compute_energy();

	check_overlaps();

	if (this->_have_us) {
		this->_update_ops();
		check_ops();
	}
	
	if (this->_overlap == true)
		this->_IO->die("There is an overlap in the initial configuration. Dying badly");
}

template<typename number>
void VMMC_CPUBackend<number>::get_settings(input_file & inp) {
	MC_CPUBackend<number>::get_settings(inp);
	int is_us, tmpi;

	if (getInputInt(&inp, "maxclust", &tmpi, 0) == KEY_FOUND) {
		_maxclust = tmpi;
		this->_IO->log(this->_IO->LOG_INFO, "Using maxclust = %i", _maxclust);
	}
	
	if (getInputInt(&inp, "small_system", &tmpi, 0) == KEY_FOUND) {
		if (tmpi > 0) {
			_small_system = true;
			this->_IO->log(this->_IO->LOG_INFO, "Using algorithm N^2 for small system");
		}
		else {
			this->_IO->log(this->_IO->LOG_INFO, "Using standard algorithm");
		}
	}
	
	if (getInputInt(&inp, "preserve_topology", &tmpi, 0) == KEY_FOUND) {
		if (tmpi > 0) {
			_preserve_topology = true;
		}
		else {
			_preserve_topology = false;
		}
	}

	if (getInputInt(&inp, "umbrella_sampling", &is_us, 0) != KEY_NOT_FOUND) {
		if (is_us > 0) {
			_have_us = true;
			getInputString(&inp, "op_file", _op_file, 1);
			getInputString(&inp, "weights_file", _weights_file, 1);
			if (getInputString(&inp, "last_hist_file", _last_hist_file, 0)
					== KEY_NOT_FOUND) {
				sprintf(_last_hist_file, "last_hist.dat");
				this->_IO->log(this->_IO->LOG_INFO,
						"Using default hist file %s", _last_hist_file);
			}
			if (getInputString(&inp, "traj_hist_file", _traj_hist_file, 0)
					== KEY_NOT_FOUND) {
				sprintf(_traj_hist_file, "traj_hist.dat");
				this->_IO->log(this->_IO->LOG_INFO,
						"Using default traj hist file %s", _traj_hist_file);
			}

			// should we reload histograms?
			if (getInputString(&inp, "init_hist_file", _init_hist_file, 0)
					== KEY_FOUND) {
				this->_IO->log(this->_IO->LOG_INFO,
						"Reloading histogram from %s", _init_hist_file);
				this->_reload_hist = true;
			} else {
				this->_reload_hist = false;
				FILE * tmpfile = fopen(_traj_hist_file, "w");
				fclose(tmpfile);
			}

			// should we extrapolate the histogram at different
			// temperatures?
			char tstring[512];
			if (getInputString(&inp, "extrapolate_hist", tstring, 0)
					== KEY_FOUND) {
				this->_IO->log(this->_IO->LOG_INFO,
						"Extrapolating temperatures .... %s", tstring);
				char * aux, deg;
				int c = 0, check;
				double * tmpt;
				tmpt = new double[100];
				aux = strtok(tstring, ",");
				while (aux != NULL) {
					//printf ("parsing %s\n", aux);
					check = sscanf(aux, "%lf %c", &(tmpt[c]), &deg);
					if (check < 1) {
						this->_IO->die(
								"Unrecognizable line in extrapolate_hist");
					}
					if (check == 1) {
						; // do nothing
					}
					if (check == 2) {
						deg = tolower(deg);
						switch (deg) {
						case 'c':
							tmpt[c] = (tmpt[c] + 273.15) * 0.1 / 300.;
							break;
						case 'k':
							tmpt[c] = tmpt[c] * 0.1 / 300.;
							break;
						default:
							this->_IO->die(
									"Unrecognizable temperature '%s' in extrapolate_hist",
									tmpt[c]);
							break;
						}
					}
					c++;
					aux = strtok(NULL, ",");
				}
				if (c == 0) {
					this->_IO->die("Nothing found in extrapolate_hist");
				} else {
					fprintf(stderr, "Extrapolating to temperatures ");
					for (int i = 0; i < c; i++)
						fprintf(stderr, "%lf ", tmpt[i]);
					fprintf(stderr, "\n");
				}
				_netemps = c;
				if (_netemps > 0) {
					_etemps = new double[_netemps];
					memcpy(_etemps, tmpt, _netemps * sizeof(double));
				}
				delete[] tmpt;
				//abort ();
			}
		} else {
			_have_us = false;
		}
	}
}

template<typename number>
inline number VMMC_CPUBackend<number>::_particle_particle_bonded_interaction_n5(Particle<number> *p, Particle<number> *q, number *stacking_en) {

	number energy = 0;

	// if pointer is NULL (0) we set its value to 0.;
	if (stacking_en != 0) {
		*stacking_en = (number) 0;
	}

	//if (p->n5 == q->index) {
		LR_vector<number> r = p->pos - q->pos;
		//LR_vector<number> r = p->pos.minimum_image(q->pos, this->_box_side);

		// FENE
		LR_vector<number> rback = r + p->pos_back - q->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;

		// we check wheter we ended up OUTSIDE of the FENE range
		if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
			this->_overlap = true;
			return (number) (1.e6 * this->_T);
		}
		energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

		// excluded volume

		// BASE-BASE
		LR_vector<number> rcenter = r + p->pos_base - q->pos_base;
		//energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
		energy += _excluded_volume_faster(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);

		// P-BASE vs. Q-BACK
		rcenter = r + p->pos_back - q->pos_base;
		//energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
		energy += _excluded_volume_faster(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);

		// P-BACK vs. Q-BASE
		rcenter = r + p->pos_base - q->pos_back;
		//energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
		energy += _excluded_volume_faster(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);

		// STACKING

		// here we define rbackghost, which is the position the backbone would be in if we were using the no-grooving model. We do this because some of the force calculations rely on this definition of rback; it's just easier than recalculating the constants in that calculation
		// NB if major-minor grooving is not in use, rback = rbackghost and everything works as it should
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> rbackghost = r + a1 * POS_BACK - b1 * POS_BACK;
		number rbackghostmod = rbackghost.module();

		LR_vector<number> rstack = r + p->pos_stack - q->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;
		
		/*
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;
		*/

		number cost4 = p->orientationT.v3 * q->orientationT.v3;
		number cost5 = p->orientationT.v3 * rstackdir;
		number cost6 =-q->orientationT.v3 * rstackdir;
		number cosphi1 = p->orientationT.v2 * rbackghost / rbackghostmod;
		number cosphi2 = q->orientationT.v2 * rbackghost / rbackghostmod;

		number f1 = this->_interaction.f1(rstackmod, STCK_F1, p->type, q->type);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5 = this->_interaction.query_mesh(-cost5,
				this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6 = this->_interaction.query_mesh(cost6,
				this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		//printf ("p->n5 = q: %d %d (%lf %lf %lf) ", p->index, q->index, rback[0], rback[1], rback[2]);
		//printf (" - %lf %lf %lf %lf %lf %lf\n", f1, f4t4, f4t5, f4t6, f5phi1, f5phi2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		if (stacking_en != 0) {
			*stacking_en = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		}
		return energy;
//
//	} else {
//		fprintf(stderr, "Wrong neighbours .._n5! %d %d %d\n", p->index, p->n5, q->index);
//		abort();
//	}
}

template<typename number>
inline number VMMC_CPUBackend<number>::_particle_particle_bonded_interaction_n3(Particle<number> *p, Particle<number> *q, number *stacking_en) {
	if (stacking_en != 0) {
		*stacking_en = 0;
	}

//	if (p->n3 == q->index) {
		number energy = (number) 0;
		
		LR_vector<number> r = q->pos - p->pos;
		//LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

		//int id1, id2;
		//id1 = 1; id2 = 0;

		// FENE
		LR_vector<number> rback = r + q->pos_back - p->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;
		
		// we check wheter we ended up OUTSIDE of the FENE range
		if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
			this->_overlap = true;
			return (number) (1.e6 * this->_T);
		}
		energy += -FENE_EPS * 0.5 * log(1. - SQR(rbackr0) / FENE_DELTA2);

		// excluded volume

		// BASE-BASE
		LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
		//energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
		energy += _excluded_volume_faster(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);

		// P-BASE vs. Q-BACK
		rcenter = r + q->pos_back - p->pos_base;
		//energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
		energy += _excluded_volume_faster(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);

		// P-BACK vs. Q-BASE
		rcenter = r + q->pos_base - p->pos_back;
		//energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
		energy += _excluded_volume_faster(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);

		// STACKING

		// here we define rbackghost, which is the position the backbone would be in if we were using the no-grooving model. We do this because some of the force calculations rely on this definition of rback; it's just easier than recalculating the constants in that calculation
		// NB if major-minor grooving is not in use, rback = rbackghost and everything works as it should
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> rbackghost = r + b1 * POS_BACK - a1 * POS_BACK;
		number rbackghostmod = rbackghost.module();

		LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;
		
		//LR_vector<number> a2 = p->orientationT.v2;
		//LR_vector<number> b2 = q->orientationT.v2;
		//LR_vector<number> a3 = p->orientationT.v3;
		//LR_vector<number> b3 = q->orientationT.v3;
		//number cost4 = a3 * b3;
		//number cost5 = a3 * rstackdir;
		//number cost6 = -b3 * rstackdir;
		//number cosphi1 = a2 * rback / rbackmod;
		//number cosphi2 = b2 * rback / rbackmod;
		number cost4 = (p->orientationT.v3 * q->orientationT.v3);
		number cost5 = p->orientationT.v3 * rstackdir;
		number cost6 =-q->orientationT.v3 * rstackdir;
		number cosphi1 = p->orientationT.v2 * rbackghost / rbackghostmod;
		number cosphi2 = q->orientationT.v2 * rbackghost / rbackghostmod;

		// functions and their derivatives needed for energies and forces
		number f1 = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5 = this->_interaction.query_mesh(-cost5,
				this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6 = this->_interaction.query_mesh(cost6,
				this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);
		
		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		if (stacking_en != 0) {
			*stacking_en = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		}
		return energy;
//	}
//	else {
//		fprintf(stderr, "Wrong neighbours .._n3! %d %d %d\n", p->index, p->n3, q->index);
//		abort();
//	}
}

template<typename number>
inline number VMMC_CPUBackend<number>::_excluded_volume_faster(
		const LR_vector<number> &r, const number sigma, const number rstar, const number b,
		const number rc) {
	number rnorm = r.norm();
	number energy = 0;

	if (rnorm < rc * rc) {
		if (rnorm > rstar * rstar) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
		} else {
			//number lj_part = SQR(sigma / rmod) * SQR(sigma / rmod)
			//		* SQR(sigma / rmod);
			number tmp = SQR(sigma)/rnorm;
			number lj_part = tmp * tmp * tmp;
			//number lj_part = (SQR(sigma)/rnorm) * (SQR(sigma)/rnorm) * (SQR(sigma)/rnorm);
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}

	return energy;
}

template<typename number>
inline number VMMC_CPUBackend<number>::_excluded_volume(
		const LR_vector<number> &r, number sigma, number rstar, number b,
		number rc) {
	number rmod = r.module();
	number energy = 0;

	if (rmod < rc) {
		if (rmod > rstar) {
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
		} else {
			number lj_part = SQR(sigma / rmod) * SQR(sigma / rmod)
					* SQR(sigma / rmod);
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}

	return energy;
}

template<typename number>
inline number VMMC_CPUBackend<number>::_particle_particle_interaction(Particle<number> *p, Particle<number> *q, number *H_energy) {

	if(H_energy != 0) {
		*H_energy = (number) 0;
	}

	//if (p->index == 26 && q->index == 3) verbose = true;
	//if (p->index == 3 && q->index == 26) verbose = true;

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);
	
	// early ejection for particles that are in the list but are
	// actually far away...
	//if (r.module() > this->_rcut) {
	if (r.norm() > this->_sqr_rcut) {
		return 0.;
	}
	
	number energy = 0;

	//bool verbose = false;
	//bool verbose = true;

	// true if p and q are Watson-Crick pairs
	//bool is_pair = (q->type + p->type == 3);
	bool is_pair = (q->btype + p->btype == 3);

	// excluded volume intaraction
	
	// BASE-BASE
	LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
	//energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
	energy += _excluded_volume_faster(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
	//if (verbose) printf (" @@ EX1: %g ", energy);

	// P-BASE vs. Q-BACK
	rcenter = r + q->pos_back - p->pos_base;
	//energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
	energy += _excluded_volume_faster(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
	//if (verbose) printf (" %g ", energy);

	// P-BACK vs. Q-BASE
	rcenter = r + q->pos_base - p->pos_back;
	//energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
	energy += _excluded_volume_faster(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
	//if (verbose) printf (" %g ", energy);

	// BACK-BACK
	rcenter = r + q->pos_back - p->pos_back;
	//energy += _excluded_volume(rcenter, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);
	energy += _excluded_volume_faster(rcenter, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);
	//if (verbose) printf (" %g\n", energy);

	// HYDROGEN BONDING
	number hb_energy = (number) 0;
	// vector, versor and magnitude of the base-base separation
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();
	LR_vector<number> rhydrodir = rhydro / rhydromod;
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		/*
		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		//LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		//LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		*/

		// angles involved in the HB interaction
		/*
		number t1 = LRACOS (-a1 * b1);
		number t2 = LRACOS (-b1 * rhydrodir);
		number t3 = LRACOS ( a1 * rhydrodir);
		number t4 = LRACOS ( a3 * b3);
		number t7 = LRACOS (-b3 * rhydrodir);
		number t8 = LRACOS ( a3 * rhydrodir);
		 */

		/*
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 =  a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 =  a3 * rhydrodir;
		*/
		number cost1 = -p->orientationT.v1 * q->orientationT.v1;
		number cost2 = -q->orientationT.v1 * rhydrodir;
		number cost3 =  p->orientationT.v1 * rhydrodir;

		number cost4 =  p->orientationT.v3 * q->orientationT.v3;
		number cost7 = -q->orientationT.v3 * rhydrodir;
		number cost8 =  p->orientationT.v3 * rhydrodir;

		// functions called at their relevant arguments
		number f1 = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
		/*
		 number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
		 number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
		 number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[HYDR_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh(cost2,
				this->_interaction.mesh_f4[HYDR_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh(cost3,
				this->_interaction.mesh_f4[HYDR_F4_THETA3]);

		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[HYDR_F4_THETA4]);
		number f4t7 = this->_interaction.query_mesh(cost7,
				this->_interaction.mesh_f4[HYDR_F4_THETA7]);
		number f4t8 = this->_interaction.query_mesh(cost8,
				this->_interaction.mesh_f4[HYDR_F4_THETA8]);
		/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
		 number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
		 number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		//if (verbose) printf ("%g %g %g %g %g %g %g @@@\n",
		//		f1, f4t1, f4t2, f4t3, f4t4, f4t7, f4t8);
		energy += hb_energy;
		//this->_U_hydr += hb_energy;
		if(H_energy != 0) {
			*H_energy = hb_energy;
		}
	}
	// END OF HYDROGEN BONDING
	
	// CROSS STACKING
	//LR_vector<number> rcstack = rhydro;
	number rcstackmod = rhydromod;
	LR_vector<number> rcstackdir = rhydrodir;
	if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
		
		// particle axes according to Allen's paper
		/*
		LR_vector<number> a1 = p->orientationT.v1;
		//LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		//LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		*/

		// angles involved in the CRST interaction
		/*
		 number t1 = LRACOS (-a1 * b1);
		 number t2 = LRACOS (-b1 * rcstackdir);
		 number t4 = LRACOS ( a3 * b3);
		 number t3 = LRACOS ( a1 * rcstackdir);
		 number t7 = LRACOS (-rcstackdir * b3);
		 number t8 = LRACOS ( rcstackdir * a3);
		 */
		/*
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 = a1 * rcstackdir;
		number cost4 = a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 = a3 * rcstackdir;
		*/
		number cost1 = -p->orientationT.v1 * q->orientationT.v1;
		number cost2 = -q->orientationT.v1 * rcstackdir;
		number cost3 =  p->orientationT.v1 * rcstackdir;
		number cost4 =  p->orientationT.v3 * q->orientationT.v3;
		number cost7 = -q->orientationT.v3 * rcstackdir;
		number cost8 =  p->orientationT.v3 * rcstackdir;

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rcstackmod, CRST_F2);
		/*
		 number f4t1 = this->_interaction.f4(t1, CRST_F4_THETA1);
		 number f4t2 = this->_interaction.f4(t2, CRST_F4_THETA2);
		 number f4t3 = this->_interaction.f4(t3, CRST_F4_THETA3);
		 number f4t4 = this->_interaction.f4(t4, CRST_F4_THETA4) + this->_interaction.f4(PI - t4, CRST_F4_THETA4);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[CRST_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh(cost2,
				this->_interaction.mesh_f4[CRST_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh(cost3,
				this->_interaction.mesh_f4[CRST_F4_THETA3]);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[CRST_F4_THETA4])
				+ this->_interaction.query_mesh(-cost4,
						this->_interaction.mesh_f4[CRST_F4_THETA4]);
		/*
		 number f4t7 = this->_interaction.f4(t7, CRST_F4_THETA7) + this->_interaction.f4(PI - t7, CRST_F4_THETA7);
		 number f4t8 = this->_interaction.f4(t8, CRST_F4_THETA8) + this->_interaction.f4(PI - t8, CRST_F4_THETA8);
		 */
		number f4t7 = this->_interaction.query_mesh(cost7,
				this->_interaction.mesh_f4[CRST_F4_THETA7])
				+ this->_interaction.query_mesh(-cost7,
						this->_interaction.mesh_f4[CRST_F4_THETA7]);
		;
		number f4t8 = this->_interaction.query_mesh(cost8,
				this->_interaction.mesh_f4[CRST_F4_THETA8])
				+ this->_interaction.query_mesh(-cost8,
						this->_interaction.mesh_f4[CRST_F4_THETA8]);
		;

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		//if (verbose) printf ("%g %g %g %g %g %g %g @@@\n",
		//		f2, f4t1, f4t2, f4t3, f4t4, f4t7, f4t8);
		energy += cstk_energy;
	}

	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
	//number rstackmod = rstack.module();
	number rstacknorm = rstack.norm();
	//if (CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
	if (SQR(CXST_RCLOW) < rstacknorm && rstacknorm < SQR(CXST_RCHIGH)) {

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;
		
		number cost4 =  a3 * b3;
		// we moved this up since it returns 0. 95% of the times in
		// a duplex, so we avoid whatever follows
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[CXST_F4_THETA4]);
		number rstackmod = sqrt(rstacknorm);
		LR_vector<number> rstackdir = rstack / rstackmod; 
			

		// angles involved in the CXST interaction
		/*
		 number t1 = LRACOS (-a1 * b1);
		 number t4 = LRACOS ( a3 * b3);
		 number t5 = LRACOS ( a3 * rstackdir);
		 number t6 = LRACOS (-b3 * rstackdir);
		 */
		number cost1 = -a1 * b1;
		//number cost4 =  a3 * b3; // moved up for eff
		number cost5 =  a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		// here we define rbackghost, which is the position the backbone would be in if we were using the no-grooving model. We do this because some of the force calculations rely on this definition of rback; it's just easier than recalculating the constants in that calculation
		// NB if major-minor grooving is not in use, rback = rbackboneghost and everything works as it should
		LR_vector<number> rbackboneghost = r + b1 * POS_BACK - a1 * POS_BACK;
		number rbackghostmod = rbackboneghost.module();
		LR_vector<number> rbackboneghostdir = rbackboneghost / rbackghostmod;
		number cosphi3 = rstackdir * (rbackboneghostdir.cross(a1));
		// old code 
		/*
		LR_vector<number> rbackbone = r + q->pos_back - p->pos_back;
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));*/
		

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rstackmod, CXST_F2);
		/*
		 number f4t1 = this->_interaction.f4(t1, CXST_F4_THETA1) + this->_interaction.f4(2 * PI - t1, CXST_F4_THETA1);
		 number f4t4 = this->_interaction.f4(t4, CXST_F4_THETA4);
		 number f4t5 = this->_interaction.f4(t5, CXST_F4_THETA5) + this->_interaction.f4(PI - t5, CXST_F4_THETA5);
		 number f4t6 = this->_interaction.f4(t6, CXST_F4_THETA6) + this->_interaction.f4(PI - t6, CXST_F4_THETA6);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[CXST_F4_THETA1]);
		//number f4t4 = this->_interaction.query_mesh(cost4,
		//		this->_interaction.mesh_f4[CXST_F4_THETA4]);
		number f4t5 = this->_interaction.query_mesh(cost5,
				this->_interaction.mesh_f4[CXST_F4_THETA5])
				+ this->_interaction.query_mesh(-cost5,
						this->_interaction.mesh_f4[CXST_F4_THETA5]);
		number f4t6 = this->_interaction.query_mesh(cost6,
				this->_interaction.mesh_f4[CXST_F4_THETA6])
				+ this->_interaction.query_mesh(-cost6,
						this->_interaction.mesh_f4[CXST_F4_THETA6]);
		number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

		number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		//printf ("%g %g %g %g %g %g @@@\n", f2, f4t1, f4t4, f4t5, f4t6, SQR(f5cosphi3));
		
		energy += cxst_energy;
	}
	
	//fflush(NULL);

	return energy;
}

inline bool find(int * clust, int size, int value) {
	int i;
	for (i = 0; i < size; i++) {
		if (clust[i] == value) {
			return true;
		}
	}
	return false;
}

/*
inline int find_i(int * clust, int size, int value) {
	int i;
	for (i = 0; i < size; i++) {
		if (clust[i] == value) {
			return i;
		}
	}
	return -1;
}
*/

template<typename number>
inline void VMMC_CPUBackend<number>::_fill_h_bonds(Particle<number> *q,
		Particle<number> *p, bool arg) {
	int n, c, tmp;

	assert (p->index != q->n3);
	assert (p->index != q->n5);
	assert (q->index != p->n3);
	assert (q->index != p->n5);

	tmp = p->get_current_neigh_index();
	p->prepare_list();
	n = p->next_neighbour();
	c = 0;
	while (n != P_VIRTUAL) {
		if (n == q->index) {
			p->h_bonds[c] = arg;
			break;
		}
		c++;
		n = p->next_neighbour();
	}
	p->set_current_neigh_index(tmp);

	tmp = q->get_current_neigh_index();
	q->prepare_list();
	n = q->next_neighbour();
	c = 0;
	while (n != P_VIRTUAL) {
		if (n == p->index) {
			q->h_bonds[c] = arg;
			break;
		}
		c++;
		n = q->next_neighbour();
	}
	q->set_current_neigh_index(tmp);

	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::_fill_e_neigh(Particle<number> *p,
		Particle<number> *q, number eij, int findex) {
	int n, c, tmp;

	assert (p->index != q->n3);
	assert (p->index != q->n5);
	assert (q->index != p->n3);
	assert (q->index != p->n5);

	p->e_neigh[findex] = eij;

	tmp = q->get_current_neigh_index();
	q->prepare_list();
	n = q->next_neighbour();
	c = 0;
	while (n != P_VIRTUAL) {
		if (n == p->index) {
			q->e_neigh[c] = eij;
			break;
		}
		c++;
		n = q->next_neighbour();
	}
	q->set_current_neigh_index(tmp);

	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::store_particle (Particle<number> * src) {
	Particle<number> *dst = &(this->_particles_old[src->index]);

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->pos_back = src->pos_back;
	dst->pos_stack = src->pos_stack;
	dst->pos_base = src->pos_base;
	dst->ext_potential = src->ext_potential;
	
	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::restore_particle (Particle<number> * dst) {
	Particle<number> *src = &(this->_particles_old[dst->index]);

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->pos_back = src->pos_back;
	dst->pos_stack = src->pos_stack;
	dst->pos_base = src->pos_base;
	dst->ext_potential = src->ext_potential;
	
	return;
}

template<typename number>
inline number VMMC_CPUBackend<number>::build_cluster_small (movestr<number> * moveptr, int maxsize, int * clust, int * size) {
	
	//printf ("before cycle...\n");	
	//if (_have_us) check_ops();
	//printf ("passed\n");

	int nclust = 1;
	clust[0] = moveptr->seed;
	Particle<number> * pp, *qq;
	number test1, test2;
	
	number pprime = 1.;

	number delta_E = 0;
	number delta_Est = 0;

	_reject_prelinks = false;
	
	set<int> prelinked_particles; //number of prelinked particles
	//set<base_pair, classcomp> poss_anomalies; // no need with N^2 algorithm
	//set<base_pair, classcomp> poss_breaks; //

	number E_anomaly = 0;
	
	number E_qq_moved;
	number E_pp_moved;
	number E_old;
	number stack_temp;
	
	// CLUSTER GENERATION
	int k = 0;
	int neigh;
	pp = &(this->_particles[clust[0]]);
	pp->inclust = true;
	
	//ppold = &(this->_particles_old[pp->index]);
	store_particle (pp);
	_move_particle(moveptr, pp);
	
	while (k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = &this->_particles[clust[k]];
		
		//printf ("recruiting from %i...\n", pp->index);

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if (_preserve_topology) {
			if (pp->pos.sqr_distance (this->_particles_old[pp->index].pos) > this->_max_move_size_sqr) {
				this->_dU = 0.;
				this->_dU_stack = 0.;
				* size = nclust;
				this->_overlap = true;
				return 0.;
			}
		}
		
		// trying to recruit bonded neighbors of pp
		if (pp->n3 != P_VIRTUAL) {
			qq = &this->_particles[pp->n3];
			if (!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3 (pp, qq, &stack_temp);
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful 
					store_particle (qq);
					_move_particle (moveptr, qq);

					E_qq_moved = _particle_particle_bonded_interaction_n3 (&(this->_particles_old[pp->index]), qq);
					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes in the cluster
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust ++;
					}
					else {
						//_r_move_particle(moveptr, qq);
						restore_particle (qq);
						prelinked_particles.insert(qq->index);
					}
				} 
				else {
					new_en3s[pp->index] = new_en5s[qq->index] = E_pp_moved;
					new_stn3s[pp->index] = new_stn5s[qq->index] = stack_temp;
					;
				}
			}
		}

		// trying to recruit 5' neighbour of pp
		if (pp->n5 != P_VIRTUAL) {
			qq = &this->_particles[pp->n5];
			if (!qq->inclust) {
				
				E_old = pp->en5;
				E_pp_moved = _particle_particle_bonded_interaction_n5 (pp, qq, &stack_temp);

				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle(moveptr, qq);
					
					E_qq_moved = _particle_particle_bonded_interaction_n5 (&this->_particles_old[pp->index], qq);

					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes to cluster
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					} 
					else {
						prelinked_particles.insert(qq->index);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
					}
				} 
				else {
					new_en5s[pp->index] = new_en3s[qq->index] = E_pp_moved;
					new_stn5s[pp->index] = new_stn3s[qq->index] = stack_temp;
					;
				}
			}
		}
		
		number tmpf = (number)0.;
		for (neigh = 0; neigh < this->_N; neigh ++) {
			qq = &this->_particles[neigh]; //qq is my neighbor
			if ((!qq->inclust) && (pp->n3 != qq->index) && (pp->n5 != qq->index) && (qq->index != pp->index)) {
				E_old = eijm_old[pp->index][qq->index];
				if (E_old == (number)0.) {
					continue;
				}
				
				E_pp_moved = _particle_particle_interaction (pp, qq, &tmpf);
				
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 >  this->_next_rand ()) {
					store_particle (qq);
					_move_particle (moveptr, qq);
					
					E_qq_moved = _particle_particle_interaction (&this->_particles_old[pp->index], qq);

					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					}
					else {
						// prelinked;
						prelinked_particles.insert(qq->index);
						//_r_move_particle (moveptr, qq);
						restore_particle (qq);
					}
				}
				else {
					//delta_E += E_pp_moved - E_old;
					eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = E_pp_moved;
					hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = (tmpf < HB_CUTOFF);
				}
			}
		}
		k ++;
	}
	*size = nclust;

	if (nclust > maxsize) {
		pprime = 0.;
		this->_dU = 0.;
		this->_dU_stack = 0.;
		this->_overlap = true;
		return pprime;
	}
	
	/*
	// Debug: print out cluster
	printf ("##@@ cluster of (%3i): ", nclust);
	for(int i = 0; i< nclust; i++) printf("%i ", clust[i]);
	printf ("\n");
	*/

	//CHECK FOR PRELINKS
	// now check if any prelinked particle is not in the cluster...
	// we reject the cluster move if we have any prelinked particles
	// that have not been fully linked at this stage
	for (int i = 0; i < nclust; i++) {
		prelinked_particles.erase(clust[i]);
	}
	if (prelinked_particles.size() > 0) {
		_reject_prelinks = true;
		this->_dU = 0;
		this->_dU_stack = 0;
		this->_overlap = true;
		return (number) 0.;
	}
	
	// fix cells...
	for (int i = 0; i < nclust; i++) {
		int old_index, new_index;
		pp = &(this->_particles[clust[i]]);
		old_index = _cells[pp->index];
		new_index = _get_cell_index(pp->pos);
		if (new_index != old_index) {
			_fix_list (pp->index, old_index, new_index);
		}
	}

	delta_E = delta_Est = E_anomaly = 0.;
	number tmpf_new, epq_new, epq_old;
	bool h_new, h_old;
	for (int i = 0; i < nclust; i++) {
		pp = &(this->_particles[clust[i]]);
		if(pp->n3 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n3]);
			if(qq->inclust == false) {
				//epq_new = _particle_particle_bonded_interaction_n3 (pp, qq, &tmpf_new);
				epq_new = new_en3s[pp->index];
				tmpf_new = new_stn3s[pp->index];
				
				delta_E += epq_new - pp->en3;
				delta_Est += tmpf_new - pp->esn3;

			}
			else {
				// in the cluster
			}
		}
		
		if(pp->n5 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n5]);
			if(qq->inclust == false) {
				//epq_new = _particle_particle_bonded_interaction_n5 (pp, qq, &tmpf_new);
				epq_new = new_en5s[pp->index];
				tmpf_new = new_stn5s[pp->index];
				
				delta_E += epq_new - pp->en5;
				delta_Est += tmpf_new - pp->esn5;
				
			}
			else {
				// in the cluster... do we care?
			}
		}

		for (neigh = 0; neigh < this->_N; neigh ++) {
			qq = &this->_particles[neigh]; //qq is my neighbor
			if ((!qq->inclust) && (pp->n3 != qq->index) && (pp->n5 != qq->index) && qq->index != pp->index) {
				epq_old = eijm_old[pp->index][qq->index];
				
				if (epq_old == (number)0.) {
					// in this case we have not considered it in the 
					// recruiting stage...
					epq_new = _particle_particle_interaction (pp, qq, &tmpf_new);
					eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = epq_new; 
					hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = (tmpf_new < HB_CUTOFF); 
				}
				else {
					// in this case, we did consider it
					epq_new = eijm[pp->index][qq->index];
				}
					
				// eijm[id1][id2] = eqp_new;
				delta_E += epq_new - epq_old;
				
				// check for anomaly of second kind;
				if (epq_old == 0. && epq_new > 0.) {
					// we have just created an overlap where there
					// was no interaction
					E_anomaly -= epq_new;
				}
				
				// check for anomaly of first kind
				if (epq_old > 0. && epq_new == 0.) {
					// we have removed an overlap
					E_anomaly += epq_old;
				}

				// fix h_bonding...
				if (_have_us) {
					h_old = hbijm_old[pp->index][qq->index];
					h_new = hbijm[pp->index][qq->index];
					//if ((pp->index == 6 && qq->index == 9) || (pp->index == 9 && qq->index == 6)) printf ("ciao.. %i %i\n", h_old, h_new);
					if (h_old != h_new) {
						if (h_old == false) {
							_op.add_hb (pp->index, qq->index);
						}
						else {
							_op.remove_hb (pp->index, qq->index);
						}
					}
				}
			}
			else { //qq in cluster
				eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = eijm_old[qq->index][pp->index];
				hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = hbijm_old[pp->index][qq->index];
			}
		}
	}

	// now we treat the case where we moved A LOT; in this case, with the
	// N^2 algorithm, no problem

	pprime *= exp((1. / this->_T) * E_anomaly);

	this->_dU = delta_E;
	this->_dU_stack = delta_Est;
	
	return pprime;
}

template<typename number>
inline number VMMC_CPUBackend<number>::build_cluster_smallish (movestr<number> * moveptr, int maxsize, int * clust, int * size) {
	
	//printf ("before cycle...\n");	
	//if (_have_us) check_ops();
	//printf ("passed\n");	

	int nclust = 1;
	clust[0] = moveptr->seed;
	Particle<number> * pp, *qq;
	number test1, test2;
	
	number pprime = 1.;

	number delta_E = 0;
	number delta_Est = 0;

	_reject_prelinks = false;
	
	set<int> prelinked_particles; //number of prelinked particles
	set<base_pair, classcomp> poss_anomalies;
	set<base_pair, classcomp> poss_breaks;

	number E_anomaly = 0;
	
	number E_qq_moved;
	number E_pp_moved;
	number E_old;
	number stack_temp;
	number H_temp;
	
	// CLUSTER GENERATION
	int k = 0;
	int icell, neigh;
	pp = &(this->_particles[clust[0]]);
	pp->inclust = true;
	
	//ppold = &(this->_particles_old[pp->index]);
	store_particle (pp);
	_move_particle(moveptr, pp);
	
	while (k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = &this->_particles[clust[k]];
		
		//printf ("recruiting from %i...\n", pp->index);

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if (_preserve_topology) {
			if (pp->pos.sqr_distance (this->_particles_old[pp->index].pos) > this->_max_move_size_sqr) {
				this->_dU = 0.;
				this->_dU_stack = 0.;
				* size = nclust;
				this->_overlap = true;
				return 0.;
			}
		}
		
		// trying to recruit bonded neighbors of pp
		if (pp->n3 != P_VIRTUAL) {
			qq = &this->_particles[pp->n3];
			if (!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3 (pp, qq, &stack_temp);
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful 
					store_particle (qq);
					_move_particle (moveptr, qq);

					E_qq_moved = _particle_particle_bonded_interaction_n3 (&(this->_particles_old[pp->index]), qq);
					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes in the cluster
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust ++;
					}
					else {
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
						prelinked_particles.insert(qq->index);
					}
				} 
				else {
					new_en3s[pp->index] = new_en5s[qq->index] = E_pp_moved;
					new_stn3s[pp->index] = new_stn5s[qq->index] = stack_temp;
					;
				}
			}
		}

		// trying to recruit 5' neighbour of pp
		if (pp->n5 != P_VIRTUAL) {
			qq = &this->_particles[pp->n5];
			if (!qq->inclust) {
				
				E_old = pp->en5;
				E_pp_moved = _particle_particle_bonded_interaction_n5 (pp, qq, &stack_temp);

				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle(moveptr, qq);
					
					E_qq_moved = _particle_particle_bonded_interaction_n5 (&this->_particles_old[pp->index], qq);

					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes to cluster
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					} 
					else {
						prelinked_particles.insert(qq->index);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
					}
				} 
				else {
					new_en5s[pp->index] = new_en3s[qq->index] = E_pp_moved;
					new_stn5s[pp->index] = new_stn3s[qq->index] = stack_temp;
					;
				}
			}
		}
		
		// a celle:
		number tmpf = (number)0.;
		for (int c = 0; c < 27; c ++) {
			icell = _neighcells[_cells[pp->index]][c];
			neigh = _vmmc_heads[icell];
			while (neigh != P_VIRTUAL) {
				qq = &this->_particles[neigh]; //qq is my neighbor
				
				if (pp->n3 == qq->index || pp->n5 == qq->index) {
					neigh = qq->_next_particle;
					continue;
				}
				
				if (qq->inclust == false) {
					E_old = eijm_old[pp->index][qq->index];
					if (E_old == (number)0.) {
						neigh = qq->_next_particle;
						continue;
					}
					
					E_pp_moved = _particle_particle_interaction (pp, qq, &tmpf);
					
					test1 = VMMC_link (E_pp_moved, E_old);
					if (test1 >  this->_next_rand ()) {
						store_particle (qq);
						_move_particle (moveptr, qq);
						
						E_qq_moved = _particle_particle_interaction (&this->_particles_old[pp->index], qq);

						test2 = VMMC_link (E_qq_moved, E_old);
						if ((test2 / test1) > this->_next_rand()) {
							clust[nclust] = qq->index;
							qq->inclust = true;
							nclust++;
						}
						else {
							// prelinked;
							prelinked_particles.insert(qq->index);
							//_r_move_particle (moveptr, qq);
							restore_particle (qq);
						}
					}
					else {
						//delta_E += E_pp_moved - E_old;
						if (E_old > (number)0) {
							poss_anomalies.insert ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						}
						if (hbijm_old[pp->index][qq->index] == true) {
							poss_breaks.insert ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						}
						eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = E_pp_moved;
						hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = (tmpf < HB_CUTOFF);
					}
				}
				neigh = qq->_next_particle;
			}
		}
		k ++;
	}
	*size = nclust;

	if (nclust > maxsize) {
		pprime = 0.;
		this->_dU = 0.;
		this->_dU_stack = 0.;
		this->_overlap = true;
		return pprime;
	}
	
	/*
	// Debug: print out cluster
	printf ("##@@ cluster of (%3i): ", nclust);
	for(int i = 0; i< nclust; i++) printf("%i ", clust[i]);
	printf ("\n");
	*/

	//CHECK FOR PRELINKS
	// now check if any prelinked particle is not in the cluster...
	// we reject the cluster move if we have any prelinked particles
	// that have not been fully linked at this stage
	for (int i = 0; i < nclust; i++) {
		prelinked_particles.erase(clust[i]);
	}
	if (prelinked_particles.size() > 0) {
		_reject_prelinks = true;
		this->_dU = 0;
		this->_dU_stack = 0;
		this->_overlap = true;
		return (number) 0.;
	}
	
	// fix cells...
	for (int i = 0; i < nclust; i++) {
		int old_index, new_index;
		pp = &(this->_particles[clust[i]]);
		old_index = _cells[pp->index];
		new_index = _get_cell_index(pp->pos);
		if (new_index != old_index) {
			_fix_list (pp->index, old_index, new_index);
		}
	}

	delta_E = delta_Est = E_anomaly = 0.;
	number tmpf_new, epq_new, epq_old;
	bool h_new, h_old;
	for (int i = 0; i < nclust; i++) {
		pp = &(this->_particles[clust[i]]);
		if(pp->n3 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n3]);
			if(qq->inclust == false) {
				//epq_new = _particle_particle_bonded_interaction_n3 (pp, qq, &tmpf_new);
				epq_new = new_en3s[pp->index];
				tmpf_new = new_stn3s[pp->index];
				
				delta_E += epq_new - pp->en3;
				delta_Est += tmpf_new - pp->esn3;

			}
			else {
				// in the cluster
			}
		}
		
		if(pp->n5 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n5]);
			if(qq->inclust == false) {
				//epq_new = _particle_particle_bonded_interaction_n5 (pp, qq, &tmpf_new);
				epq_new = new_en5s[pp->index];
				tmpf_new = new_stn5s[pp->index];
				
				delta_E += epq_new - pp->en5;
				delta_Est += tmpf_new - pp->esn5;
				
			}
			else {
				// in the cluster... do we care?
			}
		}

		for (int c = 0; c < 27; c ++) {
			icell = _neighcells[_cells[pp->index]][c];
			neigh = _vmmc_heads[icell];
			while (neigh != P_VIRTUAL) {
				qq = &(this->_particles[neigh]);

				if (pp->n3 == qq->index || pp->n5 == qq->index) {
					neigh = qq->_next_particle;
					continue;
				}
				
				if (qq->inclust == false) {
					epq_old = eijm_old[pp->index][qq->index]; // new
					h_old = hbijm_old[pp->index][qq->index];
					
					if (epq_old == (number)0.) {
						// in this case we have not considered it in the 
						// recruiting stage...
						epq_new = _particle_particle_interaction (pp, qq, &tmpf_new);
						eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = epq_new; 
						h_new = hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = (tmpf_new < HB_CUTOFF); 
					}
					else {
						// in this case, we did consider it
						epq_new = eijm[pp->index][qq->index];
						h_new = hbijm[pp->index][qq->index];
					}
						
					// eijm[id1][id2] = eqp_new;
					delta_E += epq_new - epq_old;
					
					// check for anomaly of second kind;
					if (epq_old == 0. && epq_new > 0.) {
						// we have just created an overlap where there
						// was no interaction
						E_anomaly -= epq_new;
					}
					
					// check for anomaly of first kind
					if (epq_old > 0.) {
						poss_anomalies.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						if (epq_new == 0.) {
							// we have removed an overlap
							E_anomaly += epq_old;
						}
					}

					// fix h_bonding...
					if (_have_us) {
						poss_breaks.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						if (h_old != h_new) {
							if (h_old == false) {
								_op.add_hb (pp->index, qq->index);
							}
							else {
								_op.remove_hb (pp->index, qq->index);
							}
						}
					}
				}
				else { //qq in cluster
					eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = eijm_old[qq->index][pp->index];
					hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = hbijm_old[pp->index][qq->index];
				}
				neigh = qq->_next_particle;
			}
		}
	}

	// now we treat the case where we moved A LOT; in this case, possibly a
	// hydrogen bond that was prensent between particles i (now in the
	// cluster) and j (not in the cluster) has been broken by a large move,
	// large enough that j is not in the neighborhood of i anymore. We have
	// the sets poss_breaks and poss_anomalies with these cases.
	set<base_pair>::iterator it;
	for (it = poss_breaks.begin(); it != poss_breaks.end(); it ++) {
		// if we get here, it means that the hydrogen bond between the base
		// pair has been broken, unless the particles are now both in the
		// cluster
		int i1, i2;
		i1 = (*it).first;
		i2 = (*it).second;
		pp = &this->_particles[i1];
		qq = &this->_particles[i2];
		//printf ("processing (%i, %i) at a late stage\n", pp->index, qq->index);
		if (!(pp->inclust && qq->inclust)) {
			//printf ("iao..\n");
			epq_new = _particle_particle_interaction (pp, qq, &H_temp);
			hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = (H_temp < HB_CUTOFF);
			//printf ("bong! bef: %i=%i after: %i=%i\n", hbijm_old[pp->index][qq->index], true, hbijm[qq->index][pp->index], true);
			if (H_temp > HB_CUTOFF) {
				_op.remove_hb(i1, i2);
				//printf ("removing hb between %i and %i to poss_breaks\n", pp->index, qq->index);
			} else {
				//printf ("doing nothing to (%i, %i)...\n", pp->index, qq->index);
				;
			}
		}
		else {
			//printf ("ignoring (%i, %i) since they are both in the cluster now...\n", pp->index, qq->index);
			;
		}
	}
	for (it = poss_anomalies.begin(); it != poss_anomalies.end(); it ++) {
		// if we get here, it means that the hydrogen bond between the base
		// pair has been broken, unless the particles are now both in the
		// cluster
		int i1, i2;
		i1 = (*it).first;
		i2 = (*it).second;
		pp = &this->_particles[i1];
		qq = &this->_particles[i2];
		if (!(pp->inclust && qq->inclust)) {
			epq_old = _particle_particle_interaction (&this->_particles_old[pp->index], &this->_particles_old[qq->index]);
			if (epq_old > (number)0) E_anomaly += epq_old; 
		}
	}

	pprime *= exp((1. / this->_T) * E_anomaly);

	this->_dU = delta_E;
	this->_dU_stack = delta_Est;
	
	//printf ("after cycle...\n");	
	//if (_have_us) check_ops();
	//printf ("passed\n");	
	
	return pprime;
}

template<typename number>
inline number VMMC_CPUBackend<number>::build_cluster_cells (movestr<number> * moveptr, int maxsize, int * clust, int * size) {

	int nclust = 1;
	clust[0] = moveptr->seed;
	Particle<number> * pp, *qq;
	number test1, test2;

	number pprime = 1.;

	number delta_E = 0;
	number delta_Est = 0;

	//check_ops();

	//this->_compute_energy();
	//this->_update_metainfo();
	//this->_check_metainfo();
	//this->_check_old_metainfo();

	_reject_prelinks = false;
	
	set<int> prelinked_particles; //number of prelinked particles
	set<base_pair, classcomp> poss_anomalies; //number of prelinked particles
	set<base_pair, classcomp> poss_breaks; //number of prelinked particles

	number E_anomaly = 0;
	
	number E_qq_moved;
	number E_pp_moved;
	number E_old;
	number stack_temp;
	number H_temp;
	
	// CLUSTER GENERATION
	int k = 0;
	int icell, neigh;
	pp = &(this->_particles[clust[0]]);
	pp->inclust = true;
	
	//ppold = &(this->_particles_old[pp->index]);
	store_particle (pp);
	_move_particle(moveptr, pp);
	
	while (k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = &this->_particles[clust[k]];
		
		//printf ("recruiting from %i...\n", pp->index);

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if (_preserve_topology) {
			if (pp->pos.sqr_distance (this->_particles_old[pp->index].pos) > this->_max_move_size_sqr) {
				this->_dU = 0.;
				this->_dU_stack = 0.;
				* size = nclust;
				this->_overlap = true;
				return 0.;
			}
		}
		
		// trying to recruit bonded neighbors of pp
		if (pp->n3 != P_VIRTUAL) {
			qq = &this->_particles[pp->n3];
			if (!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3 (pp, qq, &stack_temp);
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful 
					store_particle (qq);
					_move_particle(moveptr, qq);

					// in case E_pp_moved created an overlap
					this->_overlap = false;

					//_r_move_particle(moveptr, pp);
					//E_qq_moved = _particle_particle_bonded_interaction_n3(pp, qq);
					//assert (ppold->n3 == qq->index);
					E_qq_moved = _particle_particle_bonded_interaction_n3 (&(this->_particles_old[pp->index]), qq);
					//_move_particle(moveptr, pp);
					
					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes in the cluster
						
						// in case E_qq_moved created an overlap
						this->_overlap = false;
						
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust ++;
						//printf ("#FULL %i, %i\n", pp->index, qq->index);
					}
					else {
						assert (this->_overlap == false);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
						prelinked_particles.insert(qq->index);
						//printf ("#PRELINKED  %i, %i\n", pp->index, qq->index);
					}
				} 
				else {
					assert (this->_overlap == false);
					//this means qq did not get prelinked, and hence
					//interaction energy qq / pp needs to be redone; if qq
					//gets eventually added, this qq/pp energy will be set
					//to old energies
					//delta_E += E_pp_moved - E_old;
					//pp->en3 = E_pp_moved;
					//qq->en5 = E_pp_moved;
					//delta_Est += stack_temp - pp->esn3;
					//pp->esn3 = stack_temp;
					//qq->esn5 = stack_temp;
					//new_en3s[pp->index] = E_pp_moved;
					//new_stn3s[pp->index] = stack_temp;
					;
				}
			}
		}

		// trying to recruit 5' neighbour of pp
		if (pp->n5 != P_VIRTUAL) {
			qq = &this->_particles[pp->n5];
			if (!qq->inclust) {
				
				E_old = pp->en5;
				E_pp_moved = _particle_particle_bonded_interaction_n5 (pp, qq, &stack_temp);
				
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle(moveptr, qq);

					this->_overlap = false;
					
					//_r_move_particle(moveptr, pp);
					//E_qq_moved = _particle_particle_bonded_interaction_n5 (pp, qq);
					//_move_particle(moveptr, pp);
					E_qq_moved = _particle_particle_bonded_interaction_n5 (&this->_particles_old[pp->index], qq);

					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes to cluster
						this->_overlap = false;
						
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					} 
					else {
						assert (this->_overlap == false);
						prelinked_particles.insert(qq->index);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
					}
				} 
				else {
					assert (this->_overlap == false);
					//this means qq did not get prelinked, and hence
					//interaction energy qq / pp needs to be redone; if qq
					//gets eventually added, this qq/pp energy will be set
					//to old energies
					//delta_E += E_pp_moved - E_old;
					//pp->en5 = E_pp_moved;
					//qq->en3 = E_pp_moved;
					//delta_Est += stack_temp - pp->esn5;
					//pp->esn5 = stack_temp;
					//qq->esn3 = stack_temp;
					//new_en5s[pp->index] = E_pp_moved;
					//new_stn5s[pp->index] = stack_temp;
					;
				}
			}
		}
		
		// a celle:
		for (int c = 0; c < 27; c ++) {
			icell = _neighcells[_cells[pp->index]][c];
			neigh = _vmmc_heads[icell];
			while (neigh != P_VIRTUAL) {
				qq = &this->_particles[neigh]; //qq is my neighbor
				
				if (pp->n3 == qq->index || pp->n5 == qq->index) {
					neigh = qq->_next_particle;
					continue;
				}

				if (qq->inclust == false) {
					
					//E_old = en_map[MY_KEY(pp->index, qq->index)];
					//_r_move_particle (moveptr, pp);
					//E_old = _particle_particle_interaction (pp, qq, &H_temp);
					//_move_particle (moveptr, pp);
					E_old = _particle_particle_interaction (&this->_particles_old[pp->index], qq, &H_temp);
					
					if (E_old == (number)0.) {
						neigh = qq->_next_particle;
						continue;
					}
					
					E_pp_moved = _particle_particle_interaction (pp, qq);
					
					test1 = VMMC_link (E_pp_moved, E_old);
					if (test1 >  this->_next_rand ()) {
						store_particle (qq);
						_move_particle (moveptr, qq);
						
						//_r_move_particle (moveptr, pp);
						//E_qq_moved = _particle_particle_interaction (pp, qq);
						//_move_particle (moveptr, pp);
						E_qq_moved = _particle_particle_interaction (&this->_particles_old[pp->index], qq);

						test2 = VMMC_link (E_qq_moved, E_old);
						if ((test2 / test1) > this->_next_rand()) {
							clust[nclust] = qq->index;
							qq->inclust = true;
							nclust++;
						}
						else {
							// prelinked;
							prelinked_particles.insert(qq->index);
							//_r_move_particle (moveptr, qq);
							restore_particle (qq);
						}
					}
					else {
						
						//delta_E += E_pp_moved - E_old;
						
						//this is the anomalous case: positive
						//interaction energy before move, 0 energy
						//after the move
						//if(E_old > 0. && E_pp_moved == 0.) {
						//	E_anomaly += E_old; // positive if overlap removed
						//}
						// we have to check the hydrogen bonds
						// and possibly the anomalies
						if (E_old > (number)0) {
							poss_anomalies.insert ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						}
						if (H_temp < HB_CUTOFF) {
							poss_breaks.insert ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
							//printf ("adding %i and %i to poss_breaks\n", pp->index, qq->index);
						}
						//poss_breaks.insert (12);
					}
				}
				neigh = qq->_next_particle;
			}
		}
		k ++;
	}
	*size = nclust;

	if (nclust > maxsize) {
		pprime = 0.;
		this->_dU = 0.;
		this->_dU_stack = 0.;
		this->_overlap = true;
		return pprime;
	}
	
	// Debug: print out cluster
	//printf ("##@@ cluster of (%3i): ", nclust);
	//for(int i = 0; i< nclust; i++) printf("%i ", clust[i]);
	//printf ("\n");

	//CHECK FOR PRELINKS
	// now check if any prelinked particle is not in the cluster...
	// we reject the cluster move if we have any prelinked particles
	// that have not been fully linked at this stage
	for (int i = 0; i < nclust; i++) {
		prelinked_particles.erase(clust[i]);
	}
	if (prelinked_particles.size() > 0) {
		//printf ("## setting pprime = 0. because of prelinked particles..\n");
		_reject_prelinks = true;
		this->_dU = 0;
		this->_dU_stack = 0;
		this->_overlap = true;
		return (number) 0.;
	}
	
	// fix cells...
	for (int i = 0; i < nclust; i++) {
		int old_index, new_index;
		pp = &(this->_particles[clust[i]]);
		old_index = _cells[pp->index];
		new_index = _get_cell_index(pp->pos);
		if (new_index != old_index) {
			_fix_list (pp->index, old_index, new_index);
		}
	}
	
	/*	
	for (int i = 0; i < nclust; i++) {
		pp = &(this->_particles[clust[i]]);
		
		if(pp->n3 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n3]);
			if(qq->inclust) {
				if(pp->en3 != new_en3s[pp->index]) {
					//fix the total energy and the total E stacking
					//delta_E += this->_particles_old[pp->index].en3 - pp->en3;
					//delta_Est += this->_particles_old[pp->index].esn3 - pp->esn3;
					// plus old minus new
					delta_E += pp->en3 - new_en3s[pp->index];
					delta_Est += pp->esn3 - new_stn3s[pp->index];

					//new_en3s[pp->index] = 
					
					//pp->en3 = new_en3s[pp->index];
					//qq->en5 = new_en3s[pp->index];
					//pp->esn3 = new_stn3s[pp->index];
					//qq->esn5 = new_stn3s[pp->index];
				}
			}
		}
		
		//we do not check for pp->n5, because it has been takeg care of above
		//(pp is n5  neighbor of qq); We check each interacting pair only once;
		//pp-n3 is the same as qq-n5.
		
		number tmpf_old, tmpf_new, epq_new, epq_old;
		bool h_new, h_old;
		for (int c = 0; c < 27; c ++) {
			icell = _neighcells[_cells[pp->index]][c];
			neigh = _vmmc_heads[icell];
			while (neigh != P_VIRTUAL) {
				qq = &(this->_particles[neigh]);
				
				if ((pp->n3 == qq->index || pp->n5 == qq->index)) {
					neigh = qq->_next_particle;
					continue;
				}
				
				// if the particle is not in the cluster
				// we check for the formation of a positive energy
				if (qq->inclust == false) {

					epq_new = _particle_particle_interaction (pp, qq, &tmpf_new);
					if (epq_new > (number) 0.) {
						E_anomaly -= epq_new; // negative if overlap introduced
					}
					
					if (_have_us) {
						_r_move_particle (moveptr, pp);
						epq_old = _particle_particle_interaction (pp, qq, &tmpf_old);
						_move_particle (moveptr, pp);
					
						h_new = tmpf_new < HB_CUTOFF;
						h_old = tmpf_old < HB_CUTOFF;
						if (h_old != h_new) {
							if (h_old == false) {
								_op.add_hb (pp->index, qq->index);
							}
							else {
								_op.remove_hb (pp->index, qq->index);
							}
						}
					}
				}
				
				// if qq is in the cluster, but was recruited by somebody
				// else, we need to fix the new energy
				if (qq->inclust == true && pp->index < qq->index) {
					
					//if (en_map[MY_KEY(pp->index, qq->index)] > 0. && new_en_map[MY_KEY(pp->index, qq->index)] == 0.) {
					//	//it was anomalous case, we put it back
					//	E_anomaly -= en_map[MY_KEY(pp->index, qq->index)]; 
					//}
					epq_new = _particle_particle_interaction (pp, qq, &tmpf_new);
					_r_move_particle (moveptr, qq);
					_r_move_particle (moveptr, pp);
					epq_old = _particle_particle_interaction (pp, qq, &tmpf_old);
					_move_particle (moveptr, qq);
					_move_particle (moveptr, pp);
					if (epq_old > 0. && epq_new == 0.) {
						E_anomaly -= epq_old;
					}
					// remove the wrongly added term to delta_E
					//delta_E += en_map[MY_KEY(pp->index, qq->index)] - new_en_map[MY_KEY(pp->index, qq->index)];

					assert ((fabs(epq_old - epq_new) < 1.e-6));
					
					delta_E += epq_old - epq_new;
				}
				neigh = qq->_next_particle;
			} // found end of linked list for cell k
		} // end of for loop on neighbouring cells
	} // end of for loop on cluster members */

	delta_E = delta_Est = E_anomaly = 0.;
	number tmpf_old, tmpf_new, epq_new, epq_old;
	bool h_new, h_old;
	//LR_vector<number> dr, dr_old;
	for (int i = 0; i < nclust; i++) {
		pp = &(this->_particles[clust[i]]);
		if(pp->n3 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n3]);
			if(qq->inclust == false) {
				epq_new = _particle_particle_bonded_interaction_n3 (pp, qq, &tmpf_new);

				if (this->_overlap) {
					return 0.;
				}

				delta_E += epq_new - pp->en3;
				delta_Est += tmpf_new - pp->esn3;

				new_en3s[pp->index] = new_en5s[qq->index] = epq_new;
				new_stn3s[pp->index] = new_stn5s[qq->index] = tmpf_new;
			}
		}
		
		if(pp->n5 != P_VIRTUAL) {
			qq = &(this->_particles[pp->n5]);
			if(qq->inclust == false) {
				//epq_new = _particle_particle_bonded_interaction_n3 (qq, pp, &tmpf_new);
				epq_new = _particle_particle_bonded_interaction_n5 (pp, qq, &tmpf_new);

				if (this->_overlap) {
					return 0.;
				}

				delta_E += epq_new - pp->en5;
				delta_Est += tmpf_new - pp->esn5;
				
				new_en5s[pp->index] = new_en3s[qq->index] = epq_new;
				new_stn5s[pp->index] = new_stn3s[qq->index] = tmpf_new;
			}
		}

		for (int c = 0; c < 27; c ++) {
			icell = _neighcells[_cells[pp->index]][c];
			neigh = _vmmc_heads[icell];
			while (neigh != P_VIRTUAL) {
				qq = &(this->_particles[neigh]);

				if (pp->n3 == qq->index || pp->n5 == qq->index) {
					neigh = qq->_next_particle;
					continue;
				}

				if (qq->inclust == false) {
					//_r_move_particle (moveptr, pp);
					//epq_old = _particle_particle_interaction (pp, qq, &tmpf_old);
					//_move_particle (moveptr, pp);
					epq_old = _particle_particle_interaction (&this->_particles_old[pp->index], qq, &tmpf_old);
					epq_new = _particle_particle_interaction (pp, qq, &tmpf_new);
					
					delta_E += epq_new - epq_old;
					
					// check for anomaly of second kind;
					if (epq_old == 0. && epq_new > 0.) {
						// we have just created an overlap where there
						// was no interaction
						E_anomaly -= epq_new;
					}
					
					// check for anomaly of first kind
					if (epq_old > 0.) {
						if (epq_new == 0.) {
							// we have removed an overlap
							E_anomaly += epq_old;
						}
						poss_anomalies.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
					}

					// fix h_bonding...
					if (_have_us) {
						h_new = tmpf_new < HB_CUTOFF;
						h_old = tmpf_old < HB_CUTOFF;
						poss_breaks.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						//printf ("removing %i and %i from poss_breaks\n", pp->index, qq->index);
						if (h_old != h_new) {
							if (h_old == false) {
								_op.add_hb (pp->index, qq->index);
							}
							else {
								_op.remove_hb (pp->index, qq->index);
							}
						}
					}
				}
				neigh = qq->_next_particle;
			}
		}
	}

	// now we treat the case where we moved A LOT; in this case, possibly a
	// hydrogen bond that was prensent between particles i (now in the
	// cluster) and j (not in the cluster) has been broken by a large move,
	// large enough that j is not in the neighborhood of i anymore. We have
	// the sets poss_breaks and poss_anomalies with these cases.
	set<base_pair>::iterator it;
	for (it = poss_breaks.begin(); it != poss_breaks.end(); it ++) {
		// if we get here, it means that the hydrogen bond between the base
		// pair has been broken, unless the particles are now both in the
		// cluster
		int i1, i2;
		i1 = (*it).first;
		i2 = (*it).second;
		pp = &this->_particles[i1];
		qq = &this->_particles[i2];
		//printf ("processing (%i, %i) at a late stage\n", pp->index, qq->index);
		if (!(pp->inclust && qq->inclust)) {
			//printf ("iao..\n");
			epq_new = _particle_particle_interaction (pp, qq, &H_temp);
			if (H_temp > HB_CUTOFF) {
				_op.remove_hb(i1, i2);
				//printf ("removing hb between %i and %i to poss_breaks\n", pp->index, qq->index);
			} else {
				//printf ("doing nothing to (%i, %i)...\n", pp->index, qq->index);
				;
			}
		}
		else {
			//printf ("ignoring (%i, %i) since they are both in the cluster now...\n", pp->index, qq->index);
			;
		}
	}
	for (it = poss_anomalies.begin(); it != poss_anomalies.end(); it ++) {
		// if we get here, it means that the hydrogen bond between the base
		// pair has been broken, unless the particles are now both in the
		// cluster
		int i1, i2;
		i1 = (*it).first;
		i2 = (*it).second;
		pp = &this->_particles[i1];
		qq = &this->_particles[i2];
		if (!(pp->inclust && qq->inclust)) {
			//printf ("miao..\n");
			//_r_move_particle (moveptr, pp);
			//_r_move_particle (moveptr, qq);
			//epq_old = _particle_particle_interaction (pp, qq);
			//_move_particle (moveptr, pp);
			//_move_particle (moveptr, qq);
			epq_old = _particle_particle_interaction (&this->_particles_old[pp->index], &this->_particles_old[qq->index]);
			if (epq_old > (number)0) E_anomaly += epq_old; 
		}
	}

	//printf ("\n\n");

	// we fix the order parameter; it gets restored if
	// the move is rejected...
	/*
	mymap_bool::iterator it;
	for (it = new_hb_map.begin(); it != new_hb_map.end(); it ++) {
		// key: it->first, value: it->second;
		if ((*it).second != hb_map[(*it).first]) {
			int i, j;
			i = ((*it).first).first;
			j = ((*it).first).second;
			if ((*it).second == true) // if there now is HB, not before
			  _op.add_hb (i, j);
			else
			  _op.remove_hb (i, j);
		}
	}*/
	
	//int id1, id2;
	//id1 = 4;
	//id2 = 3;
//	printf ("DEBUG %i: (%i,%i) %g %g - %g %g - %g %g - % 12.10f (map % 12.10f %12.10f)\n", __LINE__,
	//printf ("DEBUG %i: (%i,%i) %g %g - %g %g - %g %g - % 12.10f (map % 12.10f)\n", __LINE__,
	//		id1, id2, (&this->_particles[id1])->en3,
	//		(&this->_particles[id2])->en5,
	//		new_en3s[id1], new_en5s[id1],
	//		new_en3s[id2], new_en5s[id2],
			//_particle_particle_bonded_interaction_n5 (
	//		_particle_particle_interaction (
	//			(&this->_particles[id1]),
	//			(&this->_particles[id2]),
	//	   	    &tmpf), 
   	//		en_map[MY_KEY(id1,id2)]);
   //			new_en_map[MY_KEY(id1,id2)]);
	//_print_pos (id1);
	//_print_pos (id2);
	//printf ("after the cicle...\n");
	//check_ops();
	//_update_ops();

	pprime *= exp((1. / this->_T) * E_anomaly);

	//assert(fabs(new_energy - original_E - delta_E) < 1.e-5 );
	this->_dU = delta_E;
	this->_dU_stack = delta_Est;
	
	return pprime;
}

/*
template<typename number>
void VMMC_CPUBackend<number>::_r_move_particle(movestr<number> * moveptr, Particle<number> *q) {
	//printf ("moving reverse %i...\n", q->index);
	if (moveptr->type == MC_MOVE_TRANSLATION) {
		q->pos -= moveptr->t;
	}
	else if (moveptr->type == MC_MOVE_ROTATION) {
		//in this case, the translation vector is the point around which we
		//rotate
		LR_vector<number> dr, drp;

		dr = q->pos.minimum_image(moveptr->t, this->_box_side);
		drp = moveptr->Rt * dr;
		q->pos += (drp - dr); // accounting for PBC
		q->orientation = moveptr->Rt * q->orientation;
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();
	}
	else {
		;
	}
	return;
}
*/
template<typename number>
inline void VMMC_CPUBackend<number>::_move_particle(movestr<number> * moveptr, Particle< number> *q) {
	//printf ("moving %i...\n", q->index);
	if (moveptr->type == MC_MOVE_TRANSLATION) {
		q->pos += moveptr->t;
	}
	else if (moveptr->type == MC_MOVE_ROTATION) {
		//in this case, the translation vector is the point around which we
		//rotate
		LR_vector<number> dr, drp;

		dr = q->pos.minimum_image(moveptr->t, this->_box_side);
		drp = moveptr->R * dr;
		q->pos += (drp - dr); // accounting for PBC
		q->orientation = moveptr->R * q->orientation;
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();
	}
	else {
		;
	}
	/*
	int old_index, new_index;
	old_index = _cells[q->index];
	new_index = _get_cell_index(q->pos);
	if (new_index != old_index) {
		_fix_list (q->index, old_index, new_index);
	}*/
	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::_fix_list (int p_index, int oldcell, int newcell) {
	int j, jold;

	// remove p_index from its old cell
	//printf ("## %i %i %i\n", oldcell, _vmmc_N_cells, _cells[p_index]);
	j = _vmmc_heads[oldcell];
	jold = P_VIRTUAL;
	assert (j != P_VIRTUAL);
	while (j != p_index) {
		jold = j;
		j = (&this->_particles[j])->_next_particle;
	}

	//printf ("j, jold: %i, %i\n", j, jold);
	//printf("no\n");
	assert (j != jold);
	assert (j != P_VIRTUAL);

	if (jold != P_VIRTUAL) {
		(&this->_particles[jold])->_next_particle = (&this->_particles[p_index])->_next_particle;
	}
	else {
		_vmmc_heads[oldcell] = (&this->_particles[p_index])->_next_particle;
	}

	// add it to the new cell
	(&this->_particles[p_index])->_next_particle = _vmmc_heads[newcell];
	_vmmc_heads[newcell] = p_index;

	_cells[p_index] = newcell;

	return;
}


template<typename number>
void VMMC_CPUBackend<number>::sim_step(llint curr_step) {
	LR_vector<number> tmp;
	//Particle<number> *pold;
	
	//printf ("checking metainfo at the beginning...\n");
	//_check_metainfo();
	//printf ("checking metainfo at the beginning: PASSED\n");

	//	_compute_energy();
	//	printf("%lf %lf\n", this->_U/this->_N, this->_U_hydr);
	//	exit(1);
	
	get_time(&this->_timer, 0);

	int * clust, nclust;
	clust = new int[this->_N];
	
	double oldweight, weight;
	int windex, oldwindex;
	oldweight = weight = 1.;
	if (_have_us)
		oldweight = _w.get_weight(_op.get_hb_states(), &oldwindex);
	
	/*
	// normalisation factor for 1/_maxclust
	number norm = log(_maxclust) + (number) 0.5772156649 + (number) (1. / 2.)
			/ (number) _maxclust + (number) (1. / 12.) / _maxclust / _maxclust;
	number logmaxclustplusone = (number) log((number) _maxclust + 1.);
	number M = (1. / norm) / (1. / (logmaxclustplusone * 2.001));
	*/

	//for (int i = 0; i < -1 ; i ++ ) {
	for (int i = 0; i < this->_N; i++) {
		if (_have_us) _op.store();
		this->_dU_stack = 0.;
		//printf ("\n##A %lf %lf \n", this->_U_stack, this->_dU_stack);
		
		// check of ext // works with TWO traps, not with one
		//number ov_c = (number) 0;
		//for (int l = 0; l < this->_N; l ++) {
		//	Particle<number> * pp = &(this->_particles[l]);
		//	pp->set_ext_potential(curr_step);
		//	ov_c += pp->ext_potential;
		//}


		// seed particle;
		int pi = (int) (drand48() * this->_N);
		Particle<number> *p = &this->_particles[pi];
		
		// this gives a random number distributed ~ 1/x (x real)
		//number trial = exp(logmaxclustplusone * drand48());
		// 1/n (n integer) is slightly different from that; to further
		// improve, we use rejection sampling (almost always accepted
		// right away).,
		/*
		number u = drand48();
		while (!(u < (1. / (norm * int(trial))) / (M * (1.
				/ (logmaxclustplusone * trial))))) {
			trial = exp(logmaxclustplusone * drand48());
			u = drand48();
		}
		maxclust = (int) trial;
		*/

		//select the move
		//printf("generating move...\n");
		movestr<number> move;
		move.seed = pi;
		move.type = (drand48() < 0.5) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;
	
		//generate translation / rotataion
		//LR_vector<number> translation;
		//LR_matrix<number> rotation;
		if (move.type == MC_MOVE_TRANSLATION) {
			move.t = LR_vector<number> (Utils::gaussian<number>(),
				Utils::gaussian<number>(), Utils::gaussian<number>()) * 
				this->_delta[MC_MOVE_TRANSLATION];
			move.R = LR_matrix<number>
				((number)1., (number) 0., (number)0.,
				 (number)0., (number) 1., (number)0.,
				 (number)0., (number) 0., (number)1.);
		}
		else {
			//translation vector is then interpreted as the axis around
			//which we rotate by move_particle() below
			//pp = &(this->_particles[clust[0]]);
			move.R = Utils::get_random_rotation_matrix_from_angle<number>
			  (this->_delta[MC_MOVE_ROTATION] * Utils::gaussian<number>());
			move.Rt = (move.R).get_transpose();
			move.t = (&this->_particles[move.seed])->pos_back + 
				(&this->_particles[move.seed])->pos;
		}
		_last_move = move.type;
		
		// build the cluster;
		//number pprime = build_cluster(pi, _maxclust, clust, &nclust, tainted, &ntainted);
		//printf("building cluster starting from %i...\n", move.seed);
		number pprime;
		if (_small_system) pprime = build_cluster_small (&move, _maxclust, clust, &nclust);
		else pprime = build_cluster_cells (&move, _maxclust, clust, &nclust);
		
		assert (nclust >=1);
		
		number delta_E_ext = 0.;
		
		// if we are not SURE to reject the move, we check the external
		// forces. Otherwise, there is no point.
		if (this->_overlap == false && pprime > 0. && nclust < this->_N) {
			//Particle<number> * q;
			for (int l = 0; l < nclust; l++) {
				p = &this->_particles[clust[l]];
				
				delta_E_ext += - p->ext_potential;
				
				p->set_ext_potential(curr_step);
				
				delta_E_ext += + p->ext_potential;

			}
			// we need to check also the tainted particles, since it is
			// possible that if we have two-body forces one was moved and
			// the other wasn't, but is in the "neighborhood" of the moved
			// particles
			pprime *= exp(-(1. / this->_T) * delta_E_ext);
		}
		
		windex = oldwindex;
		weight = oldweight;
		if (_have_us) {
			weight = _w.get_weight(&_op, &windex);
			pprime *= weight / oldweight;
		}
		
		/*
		printf ("cluster: ");
		for (int l = 0; l < nclust; l ++) {
			printf ("%i ", clust[l]);
		}
		printf ("\n");
		*/
		
		// uncomment to check the energy at a given time step.
		// may be useful for debugging purposes
		//if (curr_step > 410000 && curr_step <= 420001)
		// printf("delta_E: %lf\n", (double)delta_E);
		//printf ("### %lf\n", this->_dU_stack);
		
		this->_tries[_last_move] ++;

		//printf("## U: %lf dU: %lf, p': %lf, nclust: %d \n", this->_U, this->_dU, pprime, nclust);
		if (this->_overlap == false && pprime > drand48()) {

			assert (this->_overlap == false);

			if (nclust < _maxclust) this->_accepted[_last_move]++;
			this->_U += this->_dU;
			this->_U_stack += this->_dU_stack;

			oldweight = weight; // if (!_have_us) oldweight = weight = 1.;
			oldwindex = windex; // if (!_have_us) oldweight = weight = 1.;

			for (int l = 0; l < nclust; l ++) {
				Particle<number> * pp, * qq;
				pp = &this->_particles[clust[l]];
				if (pp->n3 != P_VIRTUAL) {
					qq = &this->_particles[pp->n3];
					if (qq->inclust == false) {
						pp->en3 = qq->en5 = new_en3s[clust[l]];
						pp->esn3 = qq->esn5 = new_stn3s[clust[l]];
					}
				}
				if (pp->n5 != P_VIRTUAL) {
					qq = &this->_particles[pp->n5];
					if (qq->inclust == false) {
						pp->en5 = qq->en3 = new_en5s[clust[l]];
						pp->esn5 = qq->esn3 = new_stn5s[clust[l]];
					}
				}

				if (_small_system) {
					// TODO: put this O(N^2)
					/*
					int icell, neigh;
					for (int c = 0; c < 27; c ++) {
						icell = _neighcells[_cells[pp->index]][c];
						neigh = _vmmc_heads[icell];
						while (neigh != P_VIRTUAL) {
							qq = &this->_particles[neigh];
							if (pp->n3 != qq->index && pp->n5 != qq->index && qq->inclust == false) {
								eijm_old[pp->index][qq->index] = eijm_old[qq->index][pp->index] = eijm[qq->index][pp->index];
								hbijm_old[pp->index][qq->index] = hbijm_old[qq->index][pp->index] = hbijm[pp->index][qq->index];
							}
							neigh = qq->_next_particle;
						}
					}*/
					for (int c = 0; c < this->_N; c ++) {
						qq = &this->_particles[c];
						if (pp->n3 != qq->index && pp->n5 != qq->index && qq->inclust == false) {
							eijm_old[pp->index][qq->index] = eijm_old[qq->index][pp->index] = eijm[qq->index][pp->index];
							hbijm_old[pp->index][qq->index] = hbijm_old[qq->index][pp->index] = hbijm[pp->index][qq->index];
						}
					}
				}
			}
			
			//printf("## accepting dU = %lf, pprime = %lf\n", this->_dU, pprime);
			//printf("## checking metainfo after accepting\n");
			//_check_metainfo();
			//if (_have_us) check_ops();
			//printf("## checking metainfo after accepting: PASSED\n");
			//_check_old_metainfo();
		}
		else {
			//move rejected
			//printf("## rejecting dU = %lf, pprime = %lf, if %i==%i just updated lists\n", this->_dU, pprime, _just_updated_lists, true);

			for (int l = 0; l < nclust; l ++) {
				int old_index, new_index;
				Particle<number> * pp;
				pp = &(this->_particles[clust[l]]);
				//_r_move_particle (&move, pp);
				restore_particle (pp);
				old_index = _cells[pp->index];
				new_index = _get_cell_index(pp->pos);
				if (new_index != old_index) {
					_fix_list (pp->index, old_index, new_index);
				}
				pp->set_ext_potential(curr_step);
			}
			
			this->_overlap = false;
			
			if (_have_us)
				_op.restore();
			//_op.print();
			
			//printf ("rejected... checking metainfo...\n");
			//if (_have_us) check_ops();
			//_check_metainfo();
			//printf ("rejected... checking metainfo: PASSED\n");
		}
		
		//check_overlaps();
		//assert (this->_overlap == false);
		
		/*	
		// check ext potential
		number c_ext = 0.;
		number c_ext_fs = 0.;
		for (int l = 0; l < this->_N; l ++) {
			Particle<number> * pp;
			pp = &(this->_particles[l]);
			c_ext += pp->ext_potential;
			pp->set_ext_potential(curr_step);
			c_ext_fs += pp->ext_potential;
		}
		if (fabs (c_ext - c_ext_fs) > 1.e-6) {
			fprintf (stderr, "%g %g -- beh\n", c_ext, c_ext_fs);
		}
		// check ext potential done*/
		
		//printf ("##G %lf %lf \n", this->_U_stack, this->_dU_stack);
		//check_energy();
		if (_have_us) {
			//assert (fabs (oldweight - _w.get_weight(&_op)) < 1.e-5);
			_h.add(oldwindex, oldweight, this->_U, this->_U_stack);
		}

		// set particles not into cluster
		for (int k = 0; k < nclust; k++) {
			(&this->_particles[clust[k]])->inclust = false;
		}
		/*
		for (int k = 0; k < ntainted; k++) {
			_tainted[tainted[k]] = false;
		}
		*/
	}

	//check_energy();
	//check_ops();
	
	delete[] clust;
	//delete[] tainted;

	get_time(&this->_timer, 1);
	process_times(&this->_timer);
}

template<typename number>
void VMMC_CPUBackend<number>::check_ops() {
	if (!_have_us) return;
	//printf ("checking OP...\n");
	assert (_have_us);
	
	int * state;
	state = (int *) malloc(_op.get_hb_parameters_count() * sizeof(int));
	memcpy(state, _op.get_hb_states(), _op.get_hb_parameters_count()
			* sizeof(int));
	_op.reset();
	
	int i, j;
	Particle<number> *p, *q;
	number hpq;
	for (i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		for (j = 0; j < i; j ++) {
			q = &this->_particles[j];
			if (p->n3 != q->index && p->n5 != q->index) {
				_particle_particle_interaction (p, q, &hpq);
				//if (p->index == 26 && q->index == 3) printf("before %g\n", epq);
				if (hpq < HB_CUTOFF) {
					//printf ("%i %i %g\n", i, j, hpq);
					_op.add_hb (p->index, q->index);
				}
			}
		}
	}

	//_print_pos(26);
	//_print_pos(3);
	//printf ("allora %lf %lf\n", _particle_particle_interaction(&(this->_particles[26]), &(this->_particles[3])), _particle_particle_interaction(&(this->_particles[26]), &(this->_particles[3]), &hpq));
	
	//_op.print ();
	int * new_state = _op.get_hb_states();
	int check = 0;
	for (i = 0; i < _op.get_hb_parameters_count(); i++) {
		if (state[i] != new_state[i])
			printf("%d should be %d \n", state[i], new_state[i]);
		check += abs(new_state[i] - state[i]);
	}

	if (check != 0) {
		printf ("CASINO\n");
		abort();
	}
	//assert (check == 0);
	//printf (" STORED\n");
	free(state);
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_update_ops() {
	//int i, j;
	
	assert (_have_us);

	_op.reset();
	
	int i, j, c;
	Particle<number> *p, *q;
	number hpq;
	for (i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		for (c = 0; c < 27; c++ ) {
			j = _vmmc_heads[_neighcells[_cells[p->index]][c]];
			while (j != P_VIRTUAL) {
				q = &this->_particles[j];
				if (j < i) {
					_particle_particle_interaction (p, q, &hpq);
					if (hpq < HB_CUTOFF) {
						_op.add_hb (i, j);
					}
				}
				j = q->_next_particle;
			}
		}
	}
	
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_check_metainfo() {
	Particle<number> *p, *q;
	number epq;
	int i, n;
	assert (this->_overlap == false);
	for (i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		n = p->n3;
		if (n != P_VIRTUAL) {
			q = &this->_particles[n];
			epq = _particle_particle_bonded_interaction_n3(p, q);
			if (this->_overlap) {
				this->_IO->die ("overlap found while checking... is your box big enough for the big rotations in your system?\n");
			}
			if (!(fabs(p->en3 - epq) < 1.e-6 || epq > 10.)) {
				printf ("LINE (part %i e %i) %i: %lf (%lf) %lf\n", p->index, q->index, __LINE__, p->en3, q->en5, epq);
				printf ("%i, %i\n", _just_updated_lists, true);
				printf ("LINE %5i; np.array([%lf, %lf, %lf]\n", __LINE__, p->pos.x, p->pos.y, p->pos.z);
				printf ("LINE %5i; np.array([%lf, %lf, %lf]\n", __LINE__, q->pos.x, q->pos.y, q->pos.z);
			}
			if (!(fabs (p->en3 - epq) < 1.e-6 || epq > 10.)) abort();
		}
		
		n = p->n5;
		if (n != P_VIRTUAL) {
			q = &this->_particles[n];
			epq = _particle_particle_bonded_interaction_n5(p, q);
			if (this->_overlap) {
				this->_IO->die ("overlap found while checking... is your box big enough for the big rotations in your system?\n");
			}
			if (!(fabs(p->en5 - epq) < 1.e-6 || epq > 10.)) {
				printf ("LINE %5i: (part %i e %i) %lf (%lf) %lf\n", __LINE__, p->index, q->index, p->en5, q->en3, epq);
				//printf ("%i, %i\n", _just_updated_lists, true);
				printf ("LINE %5i; np.array([%lf, %lf, %lf]\n", __LINE__, p->pos.x, p->pos.y, p->pos.z);
				printf ("LINE %5i; np.array([%lf, %lf, %lf]\n", __LINE__, q->pos.x, q->pos.y, q->pos.z);
			}
			if (!(fabs (p->en5 - epq) < 1.e-6 || epq > 10)) abort();
		}
		
		if (_small_system) {
			number tmpf;
			for (int j = 0; j < this->_N; j ++) {
				q = &this->_particles[j];
				if (p->n3 != q->index && p->n5 != q->index && p->index < q->index) {
					epq = _particle_particle_interaction (p, q, &tmpf);
					int id1, id2;
					id1 = MAX(p->index, q->index);
					id2 = MIN(p->index, q->index);
					if (fabs(epq - eijm_old[id1][id2]) > 1.e-5) {
						printf ("wrong value stored %i %i %g %g\n", i, j, epq, eijm_old[id1][id2]);
					}
					if (fabs(epq - eijm_old[id2][id1]) > 1.e-5) {
						printf ("wrong value stored %i %i %g %g\n", i, j, epq, eijm_old[id1][id2]);
					}
					if (_have_us) {
						if (hbijm_old[id1][id2] != (tmpf < HB_CUTOFF)) {
							printf ("wrong value stored for hb %i %i %i %g\n", i, j, hbijm_old[id1][id2], tmpf);
						}
						if (hbijm_old[id2][id1] != (tmpf < HB_CUTOFF)) {
							printf ("wrong value stored for hb %i %i %i %g\n", i, j, hbijm_old[id1][id2], tmpf);
						}
					}
					if (!(fabs (epq - eijm_old[id1][id2]) < 1.e-5 || epq > 10)) abort();
					if (!(fabs (epq - eijm_old[id2][id1]) < 1.e-5 || epq > 10)) abort();
					if (_have_us && (hbijm_old[id1][id2] != (tmpf < HB_CUTOFF))) abort();
					if (_have_us && (hbijm_old[id2][id1] != (tmpf < HB_CUTOFF))) abort();
				}
			}
		}
	}
	
	assert (this->_overlap == false);

	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_print_pos(int id) {
	Particle<number> * p;
	p = &this->_particles[id];
	printf ("%5i - np.array([% 10.6e % 10.6e % 10.6e])\n        np.array([% 10.6e % 10.6e % 10.6e])\n        np.array([% 10.6e % 10.6e % 10.6e])\n        np.array([% 10.6e % 10.6e % 10.6e])\n", id, p->pos.x, p->pos.y, p->pos.z, 
			p->orientationT.v1.x, 
			p->orientationT.v1.y, 
			p->orientationT.v1.z, 
			p->orientationT.v2.x, 
			p->orientationT.v2.y, 
			p->orientationT.v2.z,
			p->orientationT.v3.x, 
			p->orientationT.v3.y, 
			p->orientationT.v3.z);
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_update_metainfo() {
	Particle<number> *p, *q;
	number epq, tmpf;
	int i, n, c;
	for (i = 0; i < this->_N; i++) {
		p = &this->_particles[i];

		n = p->n3;
		p->en3 = p->en5 = (number) 0;
		p->esn3 = p->esn5 = (number) 0;
		if (n != P_VIRTUAL) {
			q = &this->_particles[n];
			//p->en3 = _particle_particle_interaction_pq(p, q);
			epq = _particle_particle_bonded_interaction_n3(p, q, &tmpf);
			p->en3 = epq;
			p->esn3 = tmpf;
			q->en5 = epq;
			q->esn5 = tmpf;
		}
		/*
		n = p->n5;
		if (n != P_VIRTUAL) {
			q = &this->_particles[n];
			//p->en5 = _particle_particle_interaction_pq(p, q);
			p->en5 = _particle_particle_bonded_interaction_n5(p, q);
			p->esn5 = _stacking_energy_pq(p, q);
		}
		*/
		p->prepare_list();
		n = p->next_neighbour();
		c = 0;
		while (n != P_VIRTUAL) {
			q = &this->_particles[n];
			// also updates h_bonds
			/*
			p->h_bonds[c] = (_particle_particle_hb(p, q) < HB_CUTOFF);
			p->e_neigh[c] = _particle_particle_interaction_pq(p, q);
			*/
			if (p->index < q->index) {
				epq = _particle_particle_interaction (p, q, &tmpf);
				_fill_e_neigh (p, q, epq, c);
				_fill_h_bonds (p, q, (tmpf < HB_CUTOFF));
			}
			c++;
			n = p->next_neighbour();
		}
	}
	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::check_overlaps() {
	int i, noverlaps;
	//number epq;
	Particle<number> *p, *q;
	
	noverlaps = 0;
	for (i = 0; i < this->_N; i ++) {
		p = &this->_particles[i];
		if (p->n3 != P_VIRTUAL) {
			this->_overlap = false;
			q = &this->_particles[p->n3];
			_particle_particle_bonded_interaction_n3 (p, q);
			if (this->_overlap) {
				noverlaps ++;
				LR_vector<number> rbb, r;
				r = p->pos - q->pos;
				rbb = r + p->pos_back - q->pos_back;
				printf ("### overlap %i and %i (%g, %g)\n", p->index, q->index, r.module(), rbb.module());
				this->_overlap = false;
			}
		}
	}
	assert (noverlaps == 0);
	if (noverlaps > 0) abort();
}

template<typename number>
inline void VMMC_CPUBackend<number>::check_energy() {
	//fprintf (stderr, "checking energy...\n");

	int i, j, c;
	Particle<number> * p, *q;

	// leviamoci sto dubbio
	// le interazioni NON sono simmetriche, ovvero
	// gli assert qui sotto falliscono, anche se raramente
	// e per poco. Chissa' perche'...
	/*
	for (i = 0; i < this->_N; i ++) {
		p = &this->_particles[i];
		number eij, eji;
		if (p->n3 != P_VIRTUAL) {
			q = &this->_particles[p->n3];
			eij = _particle_particle_bonded_interaction_n3 (p, q);
			eji = _particle_particle_bonded_interaction_n5 (q, p);
			if (!(fabs (eij - eji) < 1.e-5)) {
				printf ("CACCA %i %i %lf %lf\n", p->index, q->index, eij, eji);
			}
		}
		if (p->n5 != P_VIRTUAL) {
			q = &this->_particles[p->n5];
			eij = _particle_particle_bonded_interaction_n5 (p, q);
			eji = _particle_particle_bonded_interaction_n3 (q, p);
			if (!(fabs (eij - eji) < 1.e-5)) {
				printf ("CACCA %i %i %lf %lf\n", p->index, q->index, eij, eji);
			}
		}
	}*/	
	printf ("@@ %lf @@ \n", this->_U);
	
	number check = (number) 0;
	number dcheck = (number) 0;
	number estnow = (number) 0;
	for (i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		if (p->n3 != P_VIRTUAL) {
			q = &this->_particles[p->n3];
			dcheck = _particle_particle_bonded_interaction_n3(p, q, &estnow);
			check += dcheck;
			if (fabs(p->en3 - dcheck) > 1.e-6) {
				printf("%lf %lf\n", p->en3, dcheck);
				abort();
			}
			if (fabs(p->esn3 - estnow) > 1.e-6) {
				printf("%lf %lf\n", p->esn3, estnow);
				abort();
			}
		}
		if (p->n5 != P_VIRTUAL) {
			q = &this->_particles[p->n5];
			dcheck = _particle_particle_bonded_interaction_n5(p, q, &estnow);
			check += dcheck;
			if (fabs(p->en5 - dcheck) > 1.e-6) {
				printf("%lf %lf\n", p->en5, dcheck);
				abort();
			}
			if (fabs(p->esn5 - estnow) > 1.e-6) {
				printf("%lf %lf\n", p->esn5, estnow);
				abort();
			}
		}

		p->prepare_list();
		j = p->next_neighbour();
		c = 0;
		while (j != P_VIRTUAL) {
			check += p->e_neigh[c];
			bool store = p->h_bonds[c];
			if (!(fabs(_particle_particle_interaction(p, &this->_particles[j]) - p->e_neigh[c]) < 1.e-6)) {
				printf("%i %i %lf %lf\n", p->index, this->_particles[j].index, p->e_neigh[c], _particle_particle_interaction(p, &this->_particles[j]));
				int l, m;
				Particle<number> *myq;
				myq = &this->_particles[j];
				myq->prepare_list();
				l = myq->next_neighbour();
				m = 0;
				while (l != P_VIRTUAL) {
					if (l == p->index) {
						break;
					}
					m++;
					l = myq->next_neighbour();
				}
				printf("%d %d %lf %lf\n", myq->index, p->index, myq->e_neigh[m], _particle_particle_interaction(myq, p));
				abort();
			}
			if (p->h_bonds[c] != store) {
				printf ("LINE %i; maybe wrong bonds?\n", __LINE__);
				abort();
			}
			//if (p->h_bonds[c] == true) printf ("found");
			c++;
			j = p->next_neighbour();
		}
	}
	
	if (!(fabs (check / 2. - this->_U) < 1.e-5)) {
		fprintf (stderr, "%lf %lf\n", check/2., this->_U);
	}
	
	assert (fabs (check / 2. - this->_U) < 1.e-5);

	// even more checking...
	/*
	 for (i = 0; i < this->_N; i ++) {
	 p = &this->_particles[i];
	 for (j = 0; j < i; j ++) {
	 q = &this->_particles[j];
	 number eij;
	 eij = _particle_particle_hb (p, q);
	 int l, m;
	 l = find_i (p->get_verlet_list(), p->get_N_neigh (), q->index);
	 m = find_i (q->get_verlet_list(), q->get_N_neigh (), p->index);
	 if (l >= 0) assert (m >= 0);
	 assert (p->h_bonds[l] == q->h_bonds[m]);
	 if (l >= 0) {
	 assert ( (eij<(HB_CUTOFF)) == p->h_bonds[l]);
	 assert ( (eij<(HB_CUTOFF)) == q->h_bonds[m]);
	 }
	 }
	 }*/

}

template<typename number>
char * VMMC_CPUBackend<number>::get_op_state_str() {
	if (_have_us) {
		int * state = _op.get_hb_states();
		char * aux;
		aux = (char *) _state_str;
		for (int i = 0; i < _op.get_hb_parameters_count(); i++) {
			sprintf(aux, "%2d ", state[i]);
			aux = (char *) _state_str + strlen(_state_str);
		}
		sprintf(aux, "%lf", _w.get_weight(state));
		return _state_str;
	} else {
		sprintf(_state_str, " ");
		return _state_str;
	}
}

template<typename number>
void VMMC_CPUBackend<number>::print_conf(llint curr_step, bool reduced,
		bool only_last) {
	MC_CPUBackend<number>::print_conf(curr_step, reduced, only_last);
	if (_have_us) {
		if (!only_last)
			this->_h.print_to_file(_traj_hist_file, curr_step, false);
		_h.print_to_file(_last_hist_file, curr_step, true);
	}
}

template<typename number>
void VMMC_CPUBackend<number>::print_conf(llint curr_step, bool only_last) {
	MC_CPUBackend<number>::print_conf(curr_step, only_last);
	if (_have_us) {
		if (!only_last)
			this->_h.print_to_file(_traj_hist_file, curr_step, false);
		_h.print_to_file(_last_hist_file, curr_step, true);
	}
}

template<typename number>
void VMMC_CPUBackend<number>::_compute_energy () {
	
	Particle<number> * p, * q;
	this->_overlap = false;
	number res = (number) 0;
	number dres, tmpf;
	
	this->_U = this->_U_hydr = (number) 0;

	for (int i = 0; i < this->_N; i ++) {
		p = &(this->_particles[i]);
		if (p->n3 != P_VIRTUAL) {
			q = &(this->_particles[p->n3]);
			dres = _particle_particle_bonded_interaction_n3 (p, q);
			res += dres;
			this->_U += dres;
			if (this->_overlap) {
				printf ("overlap found between particle %i and %i\n", p->index, q->index);
				_print_pos (2);
				_print_pos (1);
				abort();
			}
		}
		for (int c = 0; c < 27; c ++) {
			int j = _vmmc_heads[_neighcells[_cells[p->index]][c]];
			while (j != P_VIRTUAL) {
				q = &(this->_particles[j]);
				if (p->n3 != q->index && p->n5 != q->index && p->index < q->index) {
					dres = _particle_particle_interaction (p, q, &tmpf);
					if (this->_overlap) {
						printf ("overlap found between particle %i and %i\n", p->index, q->index);
						printf ("celle: %i %i\n", _cells[p->index], _cells[q->index]);
						abort();
					}
					this->_U += dres;
					this->_U_hydr += tmpf;
				}
				j = q->_next_particle;
			}
		}
	}
	
	if (this->_overlap) {
	    this->_IO->die ("overlap found. Aborting..\n");
	}

}

template<typename number>
number VMMC_CPUBackend<number>::_compute_energy_n2 () {
	this->_overlap = false;
	Particle<number> * p, * q;
	number res = (number) 0;
	for (int i = 0; i < this->_N; i ++) {
		p = &(this->_particles[i]);
		if (p->n3 != P_VIRTUAL) {
			q = &(this->_particles[p->n3]);
			res += _particle_particle_bonded_interaction_n3 (p, q);
		}
	}

	printf ("partial... %lf\n", res);


	for (int i = 0; i < this->_N; i ++) {
		p = &(this->_particles[i]);
		for (int j = i; j < this->_N; j ++) {
			q = &(this->_particles[j]);
			if (p->n3 == q->index) continue;
			if (p->n5 == q->index) continue;
			res += _particle_particle_interaction (p, q);
		}
	}
	if (this->_overlap) printf("overlappp..\n");
	return res;
}

template<typename number>
void VMMC_CPUBackend<number>::print_energy(llint curr_step) {
	if (_have_us)
		_op.store();

	//number old_e_n2 = _compute_energy_n2();
	//number old_e = this->_U;
	//check_overlaps ();
	_check_metainfo();
	//printf ("metainfo done...\n");

	//printf ("had %g %g\n", this->_U/this->_N, this->_U_hydr);
	
	//MC_CPUBackend<number>::_compute_energy();
	_compute_energy();

	//printf("beofre printing...\n");
	if (_have_us) check_ops ();
	
	/*
	if (fabs(this->_U - old_e) > 1.e-4) {
		//printf ("LINE %5i: CACCA %lf %lf %lf\n", __LINE__, this->_U, old_e, old_e_n2);
		printf ("LINE %5i: STORED and COMPUTED ENERGIES DON'T MATCH: %lf %lf\nABORTING NOW", __LINE__, this->_U, old_e);
		abort();
	}*/

	if (_have_us)
		_op.restore();

	this->_IO->print_energy(*this, curr_step);
}

template<typename number>
void VMMC_CPUBackend<number>::_update_lists() {
	abort();
	if (this->_N > 100) {
		SimBackend<number>::_update_lists();
		return;
	}

	//printf ("eccomi->..\n");
	Particle<number> *p, *q;
	
	for(int i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		p->reset_lists();
	}

	for(int i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		for (int j = i + 1; j < this->_N; j ++) {
			q = &this->_particles[j];
			
			if (p->n3 == q->index || p->n5 == q->index) {
				continue;
			}
			
			if (p->pos.minimum_image(q->pos, this->_box_side).norm() < this->_sqr_rverlet) {
				p->add_neighbour(q->index);
				q->add_neighbour(p->index);
			}
		}
	}

	this->_are_lists_old = false;
}

template<typename number>
void VMMC_CPUBackend<number>::_init_cells() {
	
	//_vmmc_N_cells_side = (int) floor(this->_box_side / sqrt(this->_sqr_rverlet));
	_vmmc_N_cells_side = (int) floor(this->_box_side / this->_rcut + 0.1);
	//printf ("### %g %lf -> %i\n", this->_box_side, sqrt(this->_sqr_rverlet), _vmmc_N_cells_side);

	/*
	// thisi is not useful here. It's useful not to make the verlet update
	// O(27 * n_cells), which can be much more than (N^2) for a small,
	// dilute system. This makes it O(min(27 * n_cells, N^2))
	// here it's detrimental
	while(_vmmc_N_cells_side > ceil(pow(2*this->_N, 1/3.)) && _vmmc_N_cells_side > 3) {
		_vmmc_N_cells_side--;
	}
	*/
	
	if(_vmmc_N_cells_side < 3) this->_IO->die("N_cells_side (%d) must be > 2", _vmmc_N_cells_side);
	
	_vmmc_N_cells = _vmmc_N_cells_side * _vmmc_N_cells_side * _vmmc_N_cells_side;
	
	_vmmc_heads = new int[_vmmc_N_cells];
	_cells = new int[this->_N];
	_neighcells = new int * [_vmmc_N_cells];
	for (int i = 0; i < _vmmc_N_cells; i ++)
	  _neighcells[i] = new int[27];

	
	/*	
	for(int i = 0; i < _vmmc_N_cells; i++) {
		_cells[i][0] = i % _vmmc_N_cells_side;
		_cells[i][1] = i / _vmmc_N_cells_side;
		_cells[i][2] = i / (_vmmc_N_cells_side * _vmmc_N_cells_side);
		fprintf (stderr, "CELLE: %i -> (%i, %i, %i)\n", i, _cells[i][0], _cells[i][1], _cells[i][2]);
	}*/

	//fprintf (stderr, "VMMC CELL INFO: N_cells=%i, box_side=%g, N_cells_side=%i, r_cut = %g\n", _vmmc_N_cells, this->_box_side, _vmmc_N_cells_side, this->_rcut);
	
	int loop_ind[3], ind[3], nneigh;	
	for(int i = 0; i < _vmmc_N_cells; i++) {
		nneigh = 0;
		ind[0] = i % _vmmc_N_cells_side;
		ind[1] = (i / _vmmc_N_cells_side) % _vmmc_N_cells_side;
		ind[2] = i / (_vmmc_N_cells_side * _vmmc_N_cells_side);
		for(int j = -1; j < 2; j++) {
			loop_ind[0] = (ind[0] + j + _vmmc_N_cells_side) % _vmmc_N_cells_side;
			for(int k = -1; k < 2; k++) {
				loop_ind[1] = (ind[1] + k + _vmmc_N_cells_side) % _vmmc_N_cells_side;
				for(int l = -1; l < 2; l++) {
					loop_ind[2] = (ind[2] + l + _vmmc_N_cells_side) % _vmmc_N_cells_side;
					int loop_index = (loop_ind[2] * _vmmc_N_cells_side + loop_ind[1]) * _vmmc_N_cells_side + loop_ind[0];
					_neighcells[i][nneigh] = loop_index;
					nneigh ++;
					//fprintf(stderr, "CELLE: %i, (%i %i %i) -> %i (%i, %i, %i)\n", i, 
					//		ind[0], ind[1], ind[2], 
					//		loop_index, loop_ind[0], loop_ind[1], loop_ind[2]);
				}
			}
		}
		assert (nneigh == 27);
	}
	//fprintf (stderr, "cells initialised...\n");
	
	for(int i = 0; i < _vmmc_N_cells; i++)
	  _vmmc_heads[i] = P_VIRTUAL;
	
	for(int i = 0; i < this->_N; i++) {
		(&this->_particles[i])->_next_particle = P_VIRTUAL;
	}

	for(int i = 0; i < this->_N; i++) {
		Particle<number> *p = &this->_particles[i];
		int cell_index = _get_cell_index(p->pos);
		int old_head = _vmmc_heads[cell_index];
		_vmmc_heads[cell_index] = i;
		_cells[i] = cell_index;
		assert (cell_index < _vmmc_N_cells);
		p->_next_particle = old_head;
	}
	//fprintf(stderr, "cells filled...\n");
	
	// check for the cells
	for (int i = 0; i < _vmmc_N_cells; i ++) {
		//fprintf (stderr, "## i %i\n", i); 
		int j = _vmmc_heads[i]; 
		//fprintf (stderr, "## j %i\n", j); 
		if (j != P_VIRTUAL) {
	//		fprintf (stderr, "cell %5i: %i ", i, j);
			j = (&this->_particles[j])->_next_particle;
			while (j != P_VIRTUAL) {
				//fprintf (stderr, "%i ", j);
				j = (&this->_particles[j])->_next_particle;
			}
	//		fprintf (stderr, "\n");
		}
	}
	//abort();
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_delete_cells() {
	delete [] _cells;
	delete [] _vmmc_heads;
	for (int i = 0; i < this->_vmmc_N_cells; i ++) delete[] _neighcells[i];
	delete [] _neighcells;
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_create_cells() {
	return;
}

template<>
inline int VMMC_CPUBackend<float>::_get_cell_index(const LR_vector<float> &pos) {
	int res = (int) ((pos.x / _box_side - floor(pos.x / _box_side)) * (1.f - FLT_EPSILON) * _N_cells_side);
	res += _vmmc_N_cells_side * ((int) ((pos.y / _box_side - floor(pos.y / _box_side)) * (1.f - FLT_EPSILON) * _N_cells_side));
	res += _vmmc_N_cells_side * _vmmc_N_cells_side *
		((int) ((pos.z / _box_side - floor(pos.z / _box_side)) * (1.f - FLT_EPSILON) * _N_cells_side));
	return res;
}

template<>
inline int VMMC_CPUBackend<double>::_get_cell_index(const LR_vector<double> &pos) {
	int res = (int) ((pos.x / _box_side - floor(pos.x / _box_side)) * (1.f - DBL_EPSILON) * _vmmc_N_cells_side);
	res += _vmmc_N_cells_side * ((int) ((pos.y / _box_side - floor(pos.y / _box_side)) * (1.f - DBL_EPSILON) * _vmmc_N_cells_side));
	res += _vmmc_N_cells_side * _vmmc_N_cells_side *
		((int) ((pos.z / _box_side - floor(pos.z / _box_side)) * (1.f - DBL_EPSILON) * _vmmc_N_cells_side));
	return res;
}

template class VMMC_CPUBackend<float> ;
template class VMMC_CPUBackend<double> ;


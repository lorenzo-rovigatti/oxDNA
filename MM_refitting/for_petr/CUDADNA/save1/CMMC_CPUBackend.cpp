/*
 * CMMC_CPUBackend.cpp
 *
 *  Created on: 26/nov/2010
 *      Author: flavio
 */


#include "CMMC_CPUBackend.h"
#include "IOManager.h"

template<typename number> CMMC_CPUBackend<number>::CMMC_CPUBackend(IOManager *IO): MC_CPUBackend<number>(IO)
{
  //_op = NULL;
  _have_us = false;
  _netemps = 0;
  _etemps = NULL;
  _maxclust = 0;
  _e_min_cutoff = (number) -1.7f;
  _e_max_cutoff = (number) -0.1f; // values good for ssDNA
}

template<typename number>
CMMC_CPUBackend<number>::~CMMC_CPUBackend() {
	if (_particles_old != NULL) delete [] _particles_old;
	if(this->_N_updates > 0) divide_given_timing(&this->_timer, 1, this->_N / (double) this->_N_updates);
	if (_netemps > 0)
	  delete [] _etemps;
}

template<typename number>
//void CMMC_CPUBackend<number>::init(ifstream &conf_input) {
void CMMC_CPUBackend<number>::init(char conf_filename[256]) {
//	MCBackend<number>::init(conf_input);
	MCBackend<number>::init(conf_filename);

	// fix maxclust if evidently wrong
	if (_maxclust < 1 || _maxclust > this->_N) {
		this->_IO->log(this->_IO->LOG_WARNING, "maxclust wrong, setting it to N");
		_maxclust = this->_N;
	}

	if (_have_us) {
		_op.init_from_file(_op_file, this->_particles, this->_N, this->_IO);
		_w.init ((const char *)_weights_file, &_op);
		if (_reload_hist) _h.init(_init_hist_file, &_op, _etemps, _netemps);
		else _h.init (&_op, _etemps, _netemps);
		_h.set_simtemp (this->_T);
	}

	_particles_old = new Particle<number>[this->_N];
	for(int i = 0; i < this->_N; i++) {
		_particles_old[i].index = i;
		_particles_old[i].type = this->_particles[i].type;
		_particles_old[i].init(this->_max_neigh);

		// this is needed for the first _compute_energy()
		this->_particles[i].set_positions();
		this->_particles[i].orientationT = this->_particles[i].orientation.get_transpose();
	}

	if(this->_delta[MC_MOVE_TRANSLATION] * sqrt(3) > this->_verlet_skin) this->_IO->die("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");

	this->_update_lists();
	this->_update_metainfo();
	this->_compute_energy();
	if (this->_have_us) this->_update_ops();
	if (this->_overlap == true) this->_IO->die("There is an overlap in the initial configuration");
}

template<typename number>
void CMMC_CPUBackend<number>::get_settings (input_file & inp) {
	MC_CPUBackend<number>::get_settings (inp);
	int is_us, tmpi;
	double tmpf;

	if(getInputInt(&inp, "maxclust", &tmpi, 0) == KEY_FOUND) {
		_maxclust = tmpi;
		this->_IO->log(this->_IO->LOG_INFO, "Using maxclust = %i", _maxclust);
	}

	if(getInputDouble(&inp, "CMMC_e_min_cutoff", &tmpf, 0) == KEY_FOUND) {
		_e_min_cutoff = (number)tmpf;
		this->_IO->log(this->_IO->LOG_INFO, "CMMC: Using min_cutoff = %lf", _e_min_cutoff);
	}

	if(getInputDouble(&inp, "CMMC_e_max_cutoff", &tmpf, 0) == KEY_FOUND) {
		_e_max_cutoff = (number)tmpf;
		this->_IO->log(this->_IO->LOG_INFO, "CMMC: Using man_cutoff = %lf", _e_max_cutoff);
	}

	if(getInputInt(&inp, "umbrella_sampling", &is_us, 0) != KEY_NOT_FOUND) {
		if (is_us > 0) {
			_have_us = true;
			getInputString(&inp, "op_file", _op_file, 1);
			getInputString(&inp, "weights_file", _weights_file, 1);
			if (getInputString(&inp, "last_hist_file", _last_hist_file, 0) == KEY_NOT_FOUND) {
				sprintf(_last_hist_file, "last_hist.dat");
				this->_IO->log(this->_IO->LOG_INFO, "Using default hist file %s", _last_hist_file);
			}
			if (getInputString(&inp, "traj_hist_file", _traj_hist_file, 0) == KEY_NOT_FOUND) {
				sprintf(_traj_hist_file, "traj_hist.dat");
				this->_IO->log(this->_IO->LOG_INFO, "Using default traj hist file %s", _traj_hist_file);
			}

			// should we reload histograms?
			if (getInputString(&inp, "init_hist_file", _init_hist_file, 0) == KEY_FOUND) {
				this->_IO->log(this->_IO->LOG_INFO, "Reloading histogram from %s", _init_hist_file);
				this->_reload_hist = true;
			}
			else {
				this->_reload_hist = false;
				FILE * tmpfile = fopen (_traj_hist_file, "w");
				fclose (tmpfile);
			}

			// should we extrapolate the histogram at different
			// temperatures?
			char tstring[512];
			if (getInputString (&inp, "extrapolate_hist", tstring, 0) == KEY_FOUND) {
				this->_IO->log(this->_IO->LOG_INFO, "Extrapolating temperatures .... %s", tstring);
				char * aux, deg;
				int c = 0, check;
				double * tmpt;
				tmpt = new double[100];
				aux = strtok (tstring, ",");
				while (aux != NULL) {
					//printf ("parsing %s\n", aux);
					check = sscanf (aux, "%lf %c", &(tmpt[c]), &deg);
					if (check < 1) {
						this->_IO->die("Unrecognizable line in extrapolate_hist");
					}
					if (check == 1) {
					   ; // do nothing
					}
					if (check == 2) {
						deg = tolower (deg);
						switch(deg) {
							case 'c':
								tmpt[c] = (tmpt[c] + 273.15) * 0.1 / 300.;
								break;
							case 'k':
								tmpt[c] = tmpt[c] * 0.1 / 300.;
								break;
							default:
								this->_IO->die("Unrecognizable temperature '%s' in extrapolate_hist", tmpt[c]);
								break;
						}
					}
					c ++;
					aux = strtok (NULL, ",");
				}
				if (c == 0) {
					this->_IO->die ("Nothing found in extrapolate_hist");
				}
				else {
					fprintf (stderr, "Extrapolating to temperatures ");
					for (int i = 0; i < c; i ++) fprintf (stderr, "%lf ", tmpt[i]);
					fprintf (stderr, "\n");
				}
				_netemps = c;
				if (_netemps > 0) {
					_etemps = new double [_netemps];
					memcpy (_etemps, tmpt, _netemps * sizeof(double));
				}
				delete [] tmpt;
				//abort ();
			}
		}
		else {
			_have_us = false;
		}
	}
}

template<typename number>
inline number CMMC_CPUBackend<number>::_particle_particle_bonded_interaction_n5(Particle<number> *p, Particle<number> *q) {
	number res;
	if (p->n5 == q->index) {
		int tmpi = p->n3;
		p->n3 = P_VIRTUAL;
		res = this->_particle_particle_bonded_interaction (p);
		p->n3 = tmpi;
		return res;
	}
	else {
		fprintf (stderr, "Wrong neighbours .._n5! %d %d %d\n", p->index, p->n5, q->index);
		abort ();
	}
}

template<typename number>
inline number CMMC_CPUBackend<number>::_particle_particle_bonded_interaction_n3(Particle<number> *p, Particle<number> *q) {
	number res;
	if (p->n3 == q->index) {
		int tmpi = p->n5;
		p->n5 = P_VIRTUAL;
		res = this->_particle_particle_bonded_interaction (p);
		p->n5 = tmpi;
		return res;
	}
	else {
		fprintf (stderr, "Wrong neighbours .._n3! %d %d %d\n", p->index, p->n5, q->index);
		abort ();
	}
}

template<typename number>
inline number CMMC_CPUBackend<number>::_excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc) {
	number rmod = r.module();
	number energy = 0;

	if(rmod < rc) {
		if(rmod > rstar) {
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
		}
		else {
			number lj_part = SQR(sigma / rmod) * SQR(sigma / rmod) * SQR(sigma / rmod);
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}

	return energy;
}

template<typename number>
inline number CMMC_CPUBackend<number>::_particle_particle_hb(Particle<number> *p, Particle<number> *q) {
	bool is_pair = (q->btype + p->btype == 3);
	number hb_energy = (number) 0;
	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
	  	LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		//LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		//LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		  /*number t1 = LRACOS (-a1 * b1);
		number t2 = LRACOS (-b1 * rhydrodir);
		number t3 = LRACOS ( a1 * rhydrodir);
		number t4 = LRACOS ( a3 * b3);
		number t7 = LRACOS (-b3 * rhydrodir);
		number t8 = LRACOS ( a3 * rhydrodir);
		*/
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

	 	 // functions called at their relevant arguments
		number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
	  	/*
		number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
	  	number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
	  	number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
		*/
		number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[HYDR_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[HYDR_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[HYDR_F4_THETA3]);

		number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[HYDR_F4_THETA4]);
		number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[HYDR_F4_THETA7]);
		number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[HYDR_F4_THETA8]);
	  	/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
	  	number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
	  	number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		//energy += hb_energy;
		//this->_U_hydr += hb_energy;
	}
	return hb_energy;
}

template<typename number>
inline number CMMC_CPUBackend<number>::_particle_particle_interaction(Particle<number> *p, Particle<number> *q) {
	// true if p and q are Watson-Crick pairs
	//bool is_pair = (q->type + p->type == 3);
	bool is_pair = (q->btype + p->btype == 3);

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	number energy = 0;

	// excluded volume

	// BASE-BASE
	LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
	energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);

	// P-BASE vs. Q-BACK
	rcenter = r + q->pos_back - p->pos_base;
	energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);

	// P-BACK vs. Q-BASE
	rcenter = r + q->pos_base - p->pos_back;
	energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);

	// BACK-BACK
	rcenter = r + q->pos_back - p->pos_back;
	energy += _excluded_volume(rcenter, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);

	// HYDROGEN BONDING
	number hb_energy = (number) 0;
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();

	// store value
	bool hb_before = false;
	if (is_pair) {
		int l = find_i (p->get_verlet_list(), p->get_N_neigh(), q->index);
		//printf ("%d %d %d %d\n", p->index, q->index, p->get_N_neigh(), l);
		assert (l >= 0 && l < p->get_N_neigh());
		assert (p->get_N_neigh() <= p->get_max_neigh());
		hb_before = p->h_bonds[l];
		//printf ("BEFA: %d\n", hb_before);
	}
	else {
		//printf ("BEFB: %d\n", hb_before);
	}

	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
	  	LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		//LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		//LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		  /*number t1 = LRACOS (-a1 * b1);
		number t2 = LRACOS (-b1 * rhydrodir);
		number t3 = LRACOS ( a1 * rhydrodir);
		number t4 = LRACOS ( a3 * b3);
		number t7 = LRACOS (-b3 * rhydrodir);
		number t8 = LRACOS ( a3 * rhydrodir);
		*/
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

	 	 // functions called at their relevant arguments
		number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
	  	/*
		number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
	  	number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
	  	number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
		*/
		number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[HYDR_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[HYDR_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[HYDR_F4_THETA3]);

		number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[HYDR_F4_THETA4]);
		number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[HYDR_F4_THETA7]);
		number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[HYDR_F4_THETA8]);
	  	/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
	  	number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
	  	number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += hb_energy;
		this->_U_hydr += hb_energy;
	}
	// fix HB;
	bool hb_after = false;
	assert (hb_after == true || hb_after == false);
	if (is_pair) {
		hb_after = (hb_energy < HB_CUTOFF);
		assert (hb_after == true || hb_after == false);
		//printf ("%d %d\n", hb_after, false);
		//printf ("%d %d\n", hb_before, true);
		if (hb_after != hb_before) {
			_fill_h_bonds (p, q, hb_energy < HB_CUTOFF);
			if (hb_after)
			  _op.add_hb (p->index, q->index);
			else
			  _op.remove_hb (p->index, q->index);
		}
	}

	// END OF HYDROGEN BONDING

	// CROSS STACKING
	LR_vector<number> rcstack = rhydro;
	//LR_vector<number> rcstack = r + q->pos_base - p->pos_base;
	number rcstackmod = rhydromod;
	//number rcstackmod = rcstack.module();
	if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
	  	LR_vector<number> rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		//LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		//LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		/*
		number t1 = LRACOS (-a1 * b1);
		number t2 = LRACOS (-b1 * rcstackdir);
		number t4 = LRACOS ( a3 * b3);
		number t3 = LRACOS ( a1 * rcstackdir);
		number t7 = LRACOS (-rcstackdir * b3);
		number t8 = LRACOS ( rcstackdir * a3);
		*/
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 =  a1 * rcstackdir;
		number cost4 =  a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 =  a3 * rcstackdir;

	 	 // functions called at their relevant arguments
		number f2   = this->_interaction.f2(rcstackmod, CRST_F2);
		/*
	  	number f4t1 = this->_interaction.f4(t1, CRST_F4_THETA1);
	  	number f4t2 = this->_interaction.f4(t2, CRST_F4_THETA2);
	  	number f4t3 = this->_interaction.f4(t3, CRST_F4_THETA3);
	  	number f4t4 = this->_interaction.f4(t4, CRST_F4_THETA4) + this->_interaction.f4(PI - t4, CRST_F4_THETA4);
		*/
	  	number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[CRST_F4_THETA1]);
	  	number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[CRST_F4_THETA2]);
	  	number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[CRST_F4_THETA3]);
	  	number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[CRST_F4_THETA4]) + this->_interaction.query_mesh (-cost4, this->_interaction.mesh_f4[CRST_F4_THETA4]);
		/*
	  	number f4t7 = this->_interaction.f4(t7, CRST_F4_THETA7) + this->_interaction.f4(PI - t7, CRST_F4_THETA7);
	  	number f4t8 = this->_interaction.f4(t8, CRST_F4_THETA8) + this->_interaction.f4(PI - t8, CRST_F4_THETA8);
		*/
	  	number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[CRST_F4_THETA7]) + this->_interaction.query_mesh (-cost7, this->_interaction.mesh_f4[CRST_F4_THETA7]);
;
	  	number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[CRST_F4_THETA8]) + this->_interaction.query_mesh (-cost8, this->_interaction.mesh_f4[CRST_F4_THETA8]);
;

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += cstk_energy;
	}

	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
	number rstackmod = rstack.module();
	if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
	  	LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		//LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		//LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		/*
		number t1 = LRACOS (-a1 * b1);
		number t4 = LRACOS ( a3 * b3);
		number t5 = LRACOS ( a3 * rstackdir);
		number t6 = LRACOS (-b3 * rstackdir);
		*/
		number cost1 = -a1 * b1;
		number cost4 =  a3 * b3;
		number cost5 =  a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		LR_vector<number> rbackbone = r + q->pos_back - p->pos_back;
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));

	 	// functions called at their relevant arguments
		number f2   = this->_interaction.f2(rstackmod, CXST_F2);
		/*
	  	number f4t1 = this->_interaction.f4(t1, CXST_F4_THETA1) + this->_interaction.f4(2 * PI - t1, CXST_F4_THETA1);
	  	number f4t4 = this->_interaction.f4(t4, CXST_F4_THETA4);
	  	number f4t5 = this->_interaction.f4(t5, CXST_F4_THETA5) + this->_interaction.f4(PI - t5, CXST_F4_THETA5);
	  	number f4t6 = this->_interaction.f4(t6, CXST_F4_THETA6) + this->_interaction.f4(PI - t6, CXST_F4_THETA6);
		*/
	  	number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[CXST_F4_THETA1]);
	  	number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[CXST_F4_THETA4]);
	  	number f4t5 = this->_interaction.query_mesh (cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]) + this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]);
	  	number f4t6 = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]) + this->_interaction.query_mesh (-cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]);
		number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

		number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy += cxst_energy;
	}

	return energy;
}

inline bool find (int * clust, int size, int value) {
	int i;
	for (i =0; i < size; i ++) {
		if (clust[i] == value) {
			return true;
		}
	}
	return false;
}

inline int find_i (int * clust, int size, int value) {
	int i;
	for (i =0; i < size; i ++) {
		if (clust[i] == value) {
			return i;
		}
	}
	return -1;
}

template<typename number>
inline void CMMC_CPUBackend<number>::_fill_h_bonds (Particle<number> *q, Particle<number> *p, bool arg) {
	int n, c, tmp;

	assert (p->index != q->n3);
	assert (p->index != q->n5);
	assert (q->index != p->n3);
	assert (q->index != p->n5);

	tmp = p->get_current_neigh_index ();
	p->prepare_list ();
	n = p->next_neighbour ();
	c = 0;
	while (n != P_VIRTUAL) {
		if (n == q->index) {
			p->h_bonds[c] = arg;
			break;
		}
		c ++;
		n = p->next_neighbour ();
	}
	p->set_current_neigh_index (tmp);

	tmp = q->get_current_neigh_index ();
	q->prepare_list ();
	n = q->next_neighbour ();
	c = 0;
	while (n != P_VIRTUAL) {
		if (n == p->index) {
			q->h_bonds[c] = arg;
			break;
		}
		c ++;
		n = q->next_neighbour ();
	}
	q->set_current_neigh_index (tmp);

	return;
}


template<typename number>
inline void CMMC_CPUBackend<number>::_fill_e_neigh (Particle<number> *p, Particle<number> *q, number eij, int findex) {
	int n, c, tmp;

	assert (p->index != q->n3);
	assert (p->index != q->n5);
	assert (q->index != p->n3);
	assert (q->index != p->n5);

	p->e_neigh[findex] = eij;

	tmp = q->get_current_neigh_index ();
	q->prepare_list ();
	n = q->next_neighbour ();
	c = 0;
	while (n != P_VIRTUAL) {
		if (n == p->index) {
			q->e_neigh[c] = eij;
			break;
		}
		c ++;
		n = q->next_neighbour ();
	}
	q->set_current_neigh_index (tmp);

	/*
	if ((p->index == 3 && q->index == 1 ) || (p->index == 3  && q->index == 1))
	  printf("@@ %d %d %lf %lf\n", p->index, q->index, p->e_neigh[c1], q->e_neigh[c2]);
	*/
	return;
}

template<typename number>
inline number CMMC_CPUBackend<number>::_stacking_energy (Particle<number> *p) {
	number energy = (number) 0;

	if(p->n3 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n3];
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = q->pos - p->pos;

		LR_vector<number> rback = r + q->pos_back - p->pos_back;
		number rbackmod = rback.module();

		// STACKING
		LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 =-b3 * rstackdir;
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		// functions and their derivatives needed for energies and forces
		number f1     = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
		number f4t4   = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5   = this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6   = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	}

	if(p->n5 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n5];
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = p->pos - q->pos;

		// FENE
		LR_vector<number> rback = r + p->pos_back - q->pos_back;
		number rbackmod = rback.module();

		// STACKING
		LR_vector<number> rstack = r + p->pos_stack - q->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 =-b3 * rstackdir;
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		number f1     = this->_interaction.f1(rstackmod, STCK_F1, p->type, q->type);
		number f4t4   = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5   = this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6   = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	}

	return energy;
}

template<typename number>
inline number CMMC_CPUBackend<number>::_stacking_energy_n3 (Particle<number> *p, Particle<number> *q) {
	number res;
	if (p->n3 == q->index) {
		int tmpi = p->n5;
		p->n5 = P_VIRTUAL;
		res = this->_stacking_energy (p);
		p->n5 = tmpi;
		return res;
	}
	else {
		fprintf (stderr, "Wrong neighbours .._n3! %d %d %d\n", p->index, p->n5, q->index);
		abort ();
	}
}

template<typename number>
inline number CMMC_CPUBackend<number>::_stacking_energy_n5 (Particle<number> *p, Particle<number> *q) {
	number res;
	if (p->n5 == q->index) {
		int tmpi = p->n3;
		p->n3 = P_VIRTUAL;
		res = this->_stacking_energy (p);
		p->n3 = tmpi;
		return res;
	}
	else {
		fprintf (stderr, "Wrong neighbours .._n5! %d %d %d\n", p->index, p->n5, q->index);
		abort ();
	}
}

template<typename number>
inline number CMMC_CPUBackend<number>::_stacking_energy_pq (Particle<number> *p, Particle<number> *q) {
	number res = (number) 0;
	if (p->n3 == q->index) {
		res += 0.5 * _stacking_energy_n3 (p, q);
		res += 0.5 * _stacking_energy_n5 (q, p);
		return res;
	}
	else if (p->n5 == q->index) {
		res += 0.5 * _stacking_energy_n5 (p, q);
		res += 0.5 * _stacking_energy_n3 (q, p);
		return res;
	}
	else {
		return (number) 0;
		//abort ();
	}
}

// only energy between cluster and environment
template<typename number>
inline number CMMC_CPUBackend<number>::_cluster_energy(int * clust, int size, bool reuse/*, number * stack_contr*/) {
	int i, neigh;
	number res = (number) 0;
	Particle<number> *p, *q;

	for (i = 0; i < size; i ++) {
		p = &this->_particles[clust[i]];
		number dres = (number)0.;
		number ddres;

		//reuse = false;
		if (reuse) {
			// we use the old information stored
			if (p->n3 != P_VIRTUAL) {
				if (!find (clust, size, p->n3)) {
					dres += p->en3;
					//this->_dU_stack -= _stacking_energy_pq (p, &this->_particles[p->n3]);
					this->_dU_stack -= p->esn3;
				}
			}
			if (p->n5 != P_VIRTUAL) {
				if (!find (clust, size, p->n5)) {
					dres += p->en5;
					//this->_dU_stack -= _stacking_energy_pq (p, &this->_particles[p->n5]);
					this->_dU_stack -= p->esn5;
				}
			}
			p->prepare_list ();
			neigh = p->next_neighbour ();
			int k = 0;
			while (neigh != P_VIRTUAL) {
				q = &this->_particles[neigh];
				if (!find (clust, size, neigh)) {
					ddres = p->e_neigh[k];
					dres += ddres;
				}
				neigh = p->next_neighbour ();
				k ++;
			}
		}
		else {
			// we compute from scratch
			if (p->n3 != P_VIRTUAL) {
				if (!find (clust, size, p->n3)) {
					q = &this->_particles[p->n3];
					ddres = _particle_particle_interaction_pq (p, q);
					dres += ddres;
					p->en3 = ddres;
					q->en5 = ddres;
					number snow = _stacking_energy_pq (p, q);
					p->esn3 = snow;
					q->esn5 = snow;
					this->_dU_stack -= snow;
				}
			}
			if (p->n5 != P_VIRTUAL) {
				if (!find (clust, size, p->n5)) {
					q = &this->_particles[p->n5];
					ddres = _particle_particle_interaction_pq (p, q);
					dres += ddres;
					p->en5 = ddres;
					q->en3 = ddres;
					number snow = _stacking_energy_pq (p, q);
					p->esn5 = snow;
					q->esn3 = snow;
					this->_dU_stack -= snow;
				}
			}
			p->prepare_list();
			neigh = p->next_neighbour ();
			int k = 0;
			while (neigh != P_VIRTUAL) {
				q = &this->_particles[neigh];
				if (!find (clust, size, neigh)) {
					ddres = _particle_particle_interaction_pq(p, q);
					//if (reuse) assert (fabs (ddres - p->e_neigh[k]) < 1.e-7);
					dres += ddres;
					this->_fill_e_neigh (p, q, ddres, k);
				}
				neigh = p->next_neighbour ();
				k ++;
			}
		}
		res += dres;
	}
	return res;
}

template<typename number>
inline number CMMC_CPUBackend<number>::_particle_particle_interaction_pq (Particle<number> *p, Particle<number> *q) {
	number res = (number) 0;
	if (p->n3 == q->index) {
		res += 0.5 * _particle_particle_bonded_interaction_n3 (p, q);
		res += 0.5 * _particle_particle_bonded_interaction_n5 (q, p);
		return res;
	}
	else if (p->n5 == q->index) {
		res += 0.5 * _particle_particle_bonded_interaction_n5 (p, q);
		res += 0.5 * _particle_particle_bonded_interaction_n3 (q, p);
		return res;
	}
	else {
		return this->_particle_particle_interaction (p, q);
	}
}

template<typename number>
inline void CMMC_CPUBackend<number>::_translate_cluster(int * clust, int size) {
	int i;
	LR_vector<number> dr;
	Particle<number> * p;
	dr.x = (drand48() - (number)0.5)*this->_delta[MC_MOVE_TRANSLATION];
	dr.y = (drand48() - (number)0.5)*this->_delta[MC_MOVE_TRANSLATION];
	dr.z = (drand48() - (number)0.5)*this->_delta[MC_MOVE_TRANSLATION];
	for (i = 0; i < size; i ++) {
		p = &this->_particles[clust[i]];
		p->pos += dr;
		if(p->pos_list.sqr_distance(p->pos) > this->_sqr_verlet_skin) {
			this->_are_lists_old = true;
		}
	}
}

template<typename number>
inline void CMMC_CPUBackend<number>::_rotate_cluster (int * clust, int size) {
	Particle<number> *p, *q;
	number t = (drand48() - (number)0.5) * this->_delta[MC_MOVE_ROTATION];
	LR_vector<number> axis = Utils::get_random_vector<number>();

	number sintheta = sin (t);
	number costheta = cos (t);
	number olcos = ((number) 1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	assert (fabs(R.determinant() - 1.) < 1.e-10);

	p = &this->_particles[clust[0]];
	LR_vector<number> dr, drp;
	for (int j = 0; j < size; j ++) {
		q = &this->_particles[clust[j]];

		dr = q->pos.minimum_image(p->pos, this->_box_side);
		drp = R * dr;
		q->pos += (drp - dr); // accounting for PBC
		//q->pos = (R * (q->pos - p->pos)) + p->pos;

		q->orientation =  R * q->orientation;
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();

		if(q->pos_list.sqr_distance(q->pos) > this->_sqr_verlet_skin) this->_are_lists_old = true;
	}
}

template<typename number>
inline void CMMC_CPUBackend<number>::build_cluster (int seed, number cutoff, int maxsize, int * clust, int * size) {
	int nclust = 1;
	clust[0] = seed;
	int k = 0;
	Particle<number> * pp, *qq;
	//while (k < nclust /*&& nclust <= maxclust*/) {
	while (k < nclust && nclust <= maxsize) {
		pp = &this->_particles[clust[k]];
		pp->inclust = true;
		if (pp->n3 != P_VIRTUAL) {
			qq = &this->_particles[pp->n3];
			//if (!find (clust, nclust, qq->index)) {
			if (!qq->inclust) {
				//if (_particle_particle_interaction_pq (pp, qq) < cutoff) {
				if (pp->en3 < cutoff) {
					clust[nclust] = qq->index;
					qq->inclust = true;
					nclust ++;
				}
			}
		}
		if (pp->n5 != P_VIRTUAL) {
			qq = &this->_particles[pp->n5];
			//if (!find (clust, nclust, qq->index)) {
			if (!qq->inclust) {
				//if (_particle_particle_interaction_pq (pp, qq) < cutoff) {
				if (pp->en5 < cutoff) {
					clust[nclust] = qq->index;
					qq->inclust = true;
					nclust ++;
				}
			}
		}
		pp->prepare_list ();
		int neigh = pp->next_neighbour ();
		int c = 0;
		while (neigh != P_VIRTUAL) {
			qq = &this->_particles[neigh];
			//if (!find (clust, nclust, qq->index)) {
			if (!qq->inclust) {
				//if (_particle_particle_interaction_pq (pp, qq) < E_CLUST_CUTOFF) {
				if (pp->e_neigh[c] < cutoff) {
					clust[nclust] = qq->index;
					qq->inclust = true;
					nclust ++;
				}
			}
			c ++;
			neigh = pp->next_neighbour ();
		}
		k ++;
	}
	* size = nclust;

	return;
}

template<typename number>
inline bool CMMC_CPUBackend<number>::check_cluster (int seed, number cutoff, int * clust, int nclust, number * ecluster) {

	bool check = true;
	number eij = (number) 0;
	Particle<number> *p, *q;
	* ecluster = (number) 0;

	for (int j = 0; j < nclust && check; j ++) {
		p = &this->_particles[clust[j]];
		if (p->n3 != P_VIRTUAL) {
			//if (!find (clust, nclust, p->n3)) {
			q = &this->_particles[p->n3];
			if (!q->inclust) {
				eij =_particle_particle_interaction_pq (p, q);
				p->en3 = eij;
				q->en5 = eij;
				if (eij < cutoff) {
					check = false;
				}
				* ecluster += eij;

				number snow = _stacking_energy_pq (p, q);
				this->_dU_stack += snow;
				p->esn3 = snow;
				q->esn5 = snow;
			}
		}
		if (p->n5 != P_VIRTUAL) {
			//if (!find (clust, nclust, p->n5)) {
			q = &this->_particles[p->n5];
			if (!q->inclust) {
				eij =_particle_particle_interaction_pq (p, q);
				p->en5 = eij;
				q->en3 = eij;
				if (eij < cutoff) {
					check = false;
				}
				* ecluster += eij;

				//this->_dU_stack += _stacking_energy_pq (p, q);
				number snow = _stacking_energy_pq (p, q);
				this->_dU_stack += snow;
				p->esn5 = snow;
				q->esn3 = snow;
			}
		}
		p->prepare_list ();
		int neigh = p->next_neighbour ();
		int c = 0;
		while (neigh != P_VIRTUAL && check) {
			q = &this->_particles[neigh];
			//if (!find (clust, nclust, neigh)) {
			if (!q->inclust) {
				eij = _particle_particle_interaction_pq (p, q);
				this->_fill_e_neigh (p, q, eij, c);
				if (p->e_neigh[c] < cutoff) {
					check = false;
				}
				* ecluster += eij;
			}
			c ++;
			neigh = p->next_neighbour();
		}
	}
	return check;
}

template<typename number>
void CMMC_CPUBackend<number>::sim_step (llint curr_step) {
	LR_vector<number> tmp;

//	_compute_energy();
//	printf("%lf %lf\n", this->_U/this->_N, this->_U_hydr);
//	exit(1);

	get_time(&this->_timer, 0);
	int * clust, nclust, maxclust, rejclust = 0;
	clust = new int[this->_N];

	double oldweight, weight;
	int windex, oldwindex;
	oldweight = weight = 1.;
	if (_have_us)
		oldweight = _w.get_weight (_op.get_hb_states(), &oldwindex);
	if (_have_us) _op.store ();

	// normalisation factor for 1/_maxclust
	number norm = log(_maxclust) + (number) 0.5772156649 + (number) (1. / 2.) / (number) _maxclust + (number) (1. / 12.) / _maxclust / _maxclust;
	number logmaxclustplusone = (number) log ((number)_maxclust + 1.);
	number M = (1./norm)/(1./(logmaxclustplusone * 2.001));

	for(int i = 0; i < this->_N; i++) {
		if (_have_us) _op.store ();
		this->_dU_stack = 0.;
		//printf ("\n##A %lf %lf \n", this->_U_stack, this->_dU_stack);

		// seed particle;
		int pi = (int) (drand48() * this->_N);
		Particle<number> *p = &this->_particles[pi];

		number E_CLUST_CUTOFF = _e_min_cutoff + (_e_max_cutoff - _e_min_cutoff) * drand48 ();

		// maximum cluster size for this move extracted from
		// a distribution ~ 1./N
		/* // rejection sampling; slower for medium and large _maxclust
		bool nc_check = true;
		while (nc_check) {
			maxclust = (int) (drand48 () * this->_maxclust) + 1; // belongs to [1, N]
			if (drand48 () < 1./(norm * maxclust)) {
				nc_check = false;
			}
		}*/

		// this gives a random number distributed ~ 1/x (x real)
		number trial = exp (logmaxclustplusone * drand48 ());
		// 1/n (n integer) is slightly different from that; to further
		// improve, we use rejection sampling (almost always accepted
		// right away).,
		number u = drand48 ();
		while (!(u < (1./(norm * int(trial)))/(M * (1./(logmaxclustplusone * trial))))) {
			trial = exp (logmaxclustplusone * drand48 ());
			u = drand48 ();
		}
		maxclust = (int) trial;

		// build the cluster;
		build_cluster (pi, E_CLUST_CUTOFF, maxclust, clust, &nclust);

		assert (nclust >=1);
		assert (maxclust >= 1);

		// we break if the cluster is too big
		this->_tries[MC_MOVE_CLUSTER_SIZE] ++;
		if (!(nclust <= maxclust)) {
			rejclust ++;
			if (_have_us) {
				_op.restore ();
				_h.add (oldwindex, oldweight, this->_U, this->_U_stack);
			}
			// set particles not to cluster
			for (int k = 0; k < nclust; k ++) {
				(&this->_particles[clust[k]])->inclust = false;
			}
			continue;
		}
		this->_accepted[MC_MOVE_CLUSTER_SIZE] ++;

		int move = (drand48() < 0.5) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;

		assert (this->_overlap == false);
		number delta_E = -_cluster_energy (clust, nclust, true);
		assert (this->_overlap == false);

		//printf ("##B %lf %lf \n", this->_U_stack, this->_dU_stack);

		number delta_E_ext = 0.;
		for (int l = 0; l < nclust; l ++) {
			p = &this->_particles[clust[l]];
			p->set_ext_potential (curr_step);
			delta_E_ext += -p->ext_potential;
		}

		if(move == MC_MOVE_TRANSLATION) {
			for (int l = 0; l < nclust; l ++) {
				p = &this->_particles[clust[l]];
				_particles_old[clust[l]].soft_copy_from (p);
			}
			_translate_cluster(clust, nclust);
			//printf ("just translated\n");
		}
		else {
			for (int l = 0; l < nclust; l ++) {
				p = &this->_particles[clust[l]];
				_particles_old[clust[l]].soft_copy_from (p);
			}
			_rotate_cluster (clust, nclust);
		}

		//printf ("##C %lf %lf \n", this->_U_stack, this->_dU_stack);

		bool just_updated_lists = false;
		get_time(&this->_timer, 2);
		if(this->_are_lists_old == true) {
			this->_update_lists();
			this->_update_metainfo();
			for (int l = 0; l < nclust; l ++) {
				p = &this->_particles[clust[l]];
				memcpy (_particles_old[clust[l]].e_neigh, p->e_neigh, (p->get_max_neigh()) * sizeof (number));
				memcpy (_particles_old[clust[l]].h_bonds, p->h_bonds, (p->get_max_neigh()) * sizeof (bool));
			}
			this->_N_updates++;
			if (_have_us) this->_update_ops();
			//printf ("just updated lists...");
			just_updated_lists = true;
		}
		get_time(&this->_timer, 3);

		//printf ("##D %lf %lf \n", this->_U_stack, this->_dU_stack);

		// check if new cluster is the same as the old one; if not,
		// detailed balance imposes to reject the move.
		bool check = true;
		number ecluster = (number) 0;
		if (!this->_overlap) {
			check = check_cluster (pi, E_CLUST_CUTOFF, clust, nclust, &ecluster);
		}

		//printf ("##E %lf %lf \n", this->_U_stack, this->_dU_stack);

		this->_tries[MC_MOVE_NEW_BOND] ++;
		this->_accepted[MC_MOVE_NEW_BOND] ++;
		if (!check) {
			this->_accepted[MC_MOVE_NEW_BOND] --;
			this->_overlap = true;
		}

		windex = oldwindex;
		weight = oldweight;
		if (!this->_overlap) {
			//delta_E += _cluster_energy (clust, nclust, false);
			// the actual energy calculation is carried out in
			// check_cluster
			//delta_E += _cluster_energy (clust, nclust, true);
			delta_E += ecluster;

			// needs to be AFTER new energy computation
			// otherwise it does not have the updated p->h_bonds
			// arrays
			if (_have_us) weight = _w.get_weight(&_op, &windex);
			for (int l = 0; l < nclust; l ++) {
				p = &this->_particles[clust[l]];
				p->set_ext_potential (curr_step);
				delta_E_ext += p->ext_potential;
			}
		}

		// uncomment to check the energy at a given time step.
		// may be useful for debugging purposes
		//if (curr_step > 410000 && curr_step <= 420001)
		// printf("delta_E: %lf\n", (double)delta_E);
		//printf ("### %lf\n", this->_dU_stack);
		this->_tries[move] ++;
		if(!this->_overlap &&
		   ( (exp(-(delta_E + delta_E_ext) / this->_T) * (weight/oldweight)) > drand48())) {
			this->_accepted[move]++;
			//this->_U += delta_E;
			//printf ("##@ %lf %lf \n", this->_U_stack, this->_dU_stack);
			this->_U_stack += this->_dU_stack;
			//printf ("### %lf %lf \n", this->_U_stack, this->_dU_stack);
			oldweight = weight; // if (!_have_us) oldweight = weight = 1.;
			oldwindex = windex; // if (!_have_us) oldweight = weight = 1.;
			//printf("## accepting\n");
		}
		else {
			if(move == MC_MOVE_TRANSLATION) {
				for (int l = 0; l < nclust; l ++) {
					p = &this->_particles[clust[l]];
					p->soft_copy_from (&_particles_old[clust[l]]);
				}
			}
			else {
				for (int l = 0; l < nclust; l ++) {
					p = &this->_particles[clust[l]];
					p->soft_copy_from (&_particles_old[clust[l]]);
					p->set_positions();
					if(p->pos_list.sqr_distance(p->pos) > this->_sqr_verlet_skin) {
						this->_are_lists_old = true;
						//printf ("just happened...\n");
					}
				}
			}

			if (this->_are_lists_old || just_updated_lists) {
				//printf ("just done....\n");
				this->_update_lists ();
				this->_update_metainfo ();
				//if (_have_us) this->_update_ops ();
			}
			else {
				// scorro i vicini dei cluster e rimetto le cose a posto
				// per quanto riguarda e_neigh...
				Particle<number> * q;
				for (int l = 0; l < nclust; l ++) {
					p = &this->_particles[clust[l]];

					// fix en3 and en5
					if (p->n3 != P_VIRTUAL) {
						//if (!find (clust, nclust, p->n3)) {
						if (!(&this->_particles[p->n3])->inclust) {
							(&this->_particles[p->n3])->en5 = p->en3;
							(&this->_particles[p->n3])->esn5 = p->esn3;
						}
					}
					if (p->n5 != P_VIRTUAL) {
						//if (!find (clust, nclust, p->n5)) {
						if (!(&this->_particles[p->n5])->inclust) {
							(&this->_particles[p->n5])->en3 = p->en5;
							(&this->_particles[p->n5])->esn3 = p->esn5;
						}
					}

					p->prepare_list ();
					int n = p->next_neighbour();
					int c = 0;
					while (n != P_VIRTUAL) {
						q = &this->_particles[n];
						//if (!find (clust, nclust, n)) {
						if (!q->inclust) {
							_fill_e_neigh (p, q, p->e_neigh[c], c);
							_fill_h_bonds (p, q, p->h_bonds[c]);
						}
						c ++;
						n = p->next_neighbour();
					}
				}
			}
			if (_have_us) _op.restore ();
			//printf("## rejecting\n");
			//p->set_ext_potential(curr_step);
		}
		this->_overlap = false;
		//printf ("##G %lf %lf \n", this->_U_stack, this->_dU_stack);
		//check_energy();
		if (_have_us) {
			//assert (fabs (oldweight - _w.get_weight(&_op)) < 1.e-5);
			_h.add (oldwindex, oldweight, this->_U, this->_U_stack);
		}
		// set particles not to cluster
		for (int k = 0; k < nclust; k ++) {
			(&this->_particles[clust[k]])->inclust = false;
		}
	}

	//check_energy();

	//check_ops();

	delete[] clust;

	get_time(&this->_timer, 1);
	process_times(&this->_timer);
}

template<typename number>
void CMMC_CPUBackend<number>::check_ops () {
	assert (_have_us);

	int * state;
	state = (int *) malloc (_op.get_hb_parameters_count () * sizeof (int));
	memcpy (state, _op.get_hb_states(), _op.get_hb_parameters_count () * sizeof (int));
	_op.reset ();
	Particle<number> *p, *q;
	int i, n, c;
	for (i = 0; i < this->_N; i ++) {
		p = &this->_particles[i];
		p->prepare_list ();
		n = p->next_neighbour ();
		c = 0;
		while (n != P_VIRTUAL) {
			q = &this->_particles[n];
			if (p->h_bonds[c] == true) {
				if (p->index < q->index) {
					_op.add_hb (p->index, q->index);
				}
			}
			c ++;
			n = p->next_neighbour ();
		}
	}

	//_op.print ();
	int * new_state = _op.get_hb_states();
	int check = 0;
	for (i = 0; i < _op.get_hb_parameters_count (); i ++) {
		if(state[i] != new_state[i])
		    printf ("%d should be %d \n", state[i], new_state[i]);
		check += abs (new_state[i] - state[i]);
	}
	assert (check == 0);
	//printf (" STORED\n");
	free (state);
	return;
}


template<typename number>
void CMMC_CPUBackend<number>::_update_ops () {
	assert (_have_us);
	_op.reset ();

	Particle<number> *p, *q;
	int i, n, c;
	for (i = 0; i < this->_N; i ++) {
		p = &this->_particles[i];
		p->prepare_list ();
		n = p->next_neighbour ();
		c = 0;
		while (n != P_VIRTUAL) {
			q = &this->_particles[n];
			if (p->h_bonds[c] == true) {
				if (p->index < q->index) {
					_op.add_hb (p->index, q->index);
				}
			}
			c ++;
			n = p->next_neighbour ();
		}
	}

	/*
	Particle<number> *p, *q;
	int i, n, c, count = 0;
	for (i = 0; i < this->_N; i ++) {
		p = &this->_particles[i];
		p->prepare_list ();
		n = p->next_neighbour ();
		c = 0;
		while (n != P_VIRTUAL) {
			q = &this->_particles[n];
			if (_particle_particle_hb (p, q) < HB_CUTOFF) {
			   printf ("found between part %d and %d\n", p->index, q->index);
			   count ++;
			}
			n = p->next_neighbour();
		}
	}
	printf ("## count: %d (%d)\n", count / 2, count);
	*/
	//_op.print ();

	//abort ();
	return;
}

template<typename number>
void CMMC_CPUBackend<number>::_update_metainfo () {
	Particle<number> *p, *q;
	int i, n, c;
	for (i = 0; i < this->_N; i ++) {
		p = &this->_particles[i];

		n = p->n3;
		p->en3 = p->en5 = (number) 0;
		p->esn3 = p->esn5 = (number) 0;
		if (n != P_VIRTUAL) {
			q = &this->_particles[n];
			p->en3 = _particle_particle_interaction_pq (p, q);
			p->esn3 = _stacking_energy_pq (p, q);
		}
		n = p->n5;
		if (n != P_VIRTUAL) {
			q = &this->_particles[n];
			p->en5 = _particle_particle_interaction_pq (p, q);
			p->esn5 = _stacking_energy_pq (p, q);
		}

		p->prepare_list ();
		n = p->next_neighbour ();
		c = 0;
		while (n != P_VIRTUAL) {
			q = &this->_particles[n];
			// also updates h_bonds
			p->h_bonds[c] = (_particle_particle_hb (p, q) < HB_CUTOFF);
			p->e_neigh[c] = _particle_particle_interaction_pq (p, q);
			c ++;
			n = p->next_neighbour ();
		}
	}
	return;
}

template<typename number>
inline void CMMC_CPUBackend<number>::check_energy() {
	//fprintf (stderr, "checking energy...\n");

	int i, j, c;
	Particle<number> * p, *q;

	/*
	// leviamoci sto dubbio
	// le interazioni NON sono simmetriche, ovvero
	// gli assert qui sotto falliscono, anche se raramente
	// e per poco. Chissa' perche'...
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

	number check = (number) 0;
	number dcheck = (number) 0;
	for (i = 0; i < this->_N; i ++)
	  {
		p = &this->_particles[i];
		if (p->n3 != P_VIRTUAL) {
			q = &this->_particles[p->n3];
			dcheck = _particle_particle_interaction_pq (p, q);
			check += dcheck;
			if (fabs (p->en3 - dcheck) > 1.e-6) {
				printf ("%lf %lf\n", p->en3, dcheck);
				abort();
			}
			if (fabs (p->esn3 - _stacking_energy_pq (p, q)) > 1.e-6) {
				printf ("%lf %lf\n", p->esn3, _stacking_energy_pq (p, q));
				abort();
			}
		}
		if (p->n5 != P_VIRTUAL) {
			q = &this->_particles[p->n5];
			dcheck = _particle_particle_interaction_pq (p, q);
			check += dcheck;
			if (fabs (p->en5 - dcheck) > 1.e-6) {
				printf ("%lf %lf\n", p->en5, dcheck);
				abort();
			}
			if (fabs (p->esn5 - _stacking_energy_pq (p, q)) > 1.e-6) {
				printf ("%lf %lf\n", p->esn5, _stacking_energy_pq (p, q));
				abort();
			}
		}

		p->prepare_list ();
		j = p->next_neighbour ();
		c = 0;
		while (j != P_VIRTUAL) {
			check += p->e_neigh[c];
			bool store = p->h_bonds[c];
			if (!(fabs (_particle_particle_interaction_pq (p, &this->_particles[j]) - p->e_neigh[c]) < 1.e-5)) {
				printf ("%i %i %lf %lf\n", p->index, this->_particles[j].index, p->e_neigh[c], _particle_particle_interaction_pq (p, &this->_particles[j]));
				int l, m;
				Particle <number> *myq;
				myq = &this->_particles[j];
				myq->prepare_list ();
				l = myq->next_neighbour();
				m = 0;
				while (l != P_VIRTUAL)
				  {
					if (l == p->index)
					  {
						break;
					  }
					m ++;
					l = myq->next_neighbour ();
				  }
				printf ("%d %d %lf %lf\n", myq->index, p->index, myq->e_neigh[m],_particle_particle_interaction_pq (myq, p));
				abort();
			}
			if (p->h_bonds[c] != store) {
				abort ();
			}
			//if (p->h_bonds[c] == true) printf ("found");
			c ++;
			j = p->next_neighbour ();
		}
	  }
	assert (fabs (check / 2. - this->_U) < 1.e-4);

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
char * CMMC_CPUBackend<number>::get_op_state_str () {
	if (_have_us) {
		int * state = _op.get_hb_states();
		char * aux;
		aux = (char *) _state_str;
		for (int i = 0; i < _op.get_hb_parameters_count(); i ++) {
			sprintf (aux, "%2d ", state[i]);
			aux = (char *) _state_str + strlen (_state_str);
		}
		sprintf (aux, "%lf", _w.get_weight(state));
		return _state_str;
	}
	else {
		sprintf (_state_str, " ");
		return _state_str;
	}
}


template<typename number>
void CMMC_CPUBackend<number>::print_conf (llint curr_step, bool reduced, bool only_last) {
	MC_CPUBackend<number>::print_conf (curr_step, reduced, only_last);
	if (_have_us) {
		if (!only_last) this->_h.print_to_file (_traj_hist_file, curr_step, false);
		_h.print_to_file (_last_hist_file, curr_step, true);
	}
}

template<typename number>
void CMMC_CPUBackend<number>::print_conf (llint curr_step, bool only_last) {
	MC_CPUBackend<number>::print_conf (curr_step, only_last);
	if (_have_us) {
		if (!only_last) this->_h.print_to_file (_traj_hist_file, curr_step, false);
		_h.print_to_file (_last_hist_file, curr_step, true);
	}
}

template<typename number>
void CMMC_CPUBackend<number>::print_energy(llint curr_step) {
	if (_have_us) _op.store();
	//_compute_energy();

	/*
		number res = (number) 0;
		number res2 = (number) 0;
	  {
		int i, j;
		for (i = 0; i < this->_N; i ++) {
			res += _stacking_energy (&this->_particles[i]);
		}
		res /= 2.;
		for (i = 0; i < this->_N; i ++) {
			for (j = 0; j < i; j ++) {
				res2 += _stacking_energy_pq (&this->_particles[i], &this->_particles[j]);
			}
		}
	  }
	number gick = this->_U;*/

	//fprintf (stderr, "@@ %lf %lf %lf %lf %lf\n", this->_U, this->_U_hydr, this->_U_stack, res, res2);
	MC_CPUBackend<number>::_compute_energy();
	/*
	if (curr_step > 10) {
	if (fabs (this->_U_stack - res2) > 1.e-6) {
	   fprintf(stderr, "here");
	   abort();
	}
	if (fabs (this->_U - gick) > 1.e-5) {
	   fprintf(stderr, "there %lf %lf", this->_U, gick);
   	   abort();
	}
	}*/
	//fprintf (stderr, "@@ %lf %lf\n", this->_U_hydr, this->_U_stack);
	if (_have_us) _op.restore();
	this->_IO->print_energy(*this, curr_step);
}

template class CMMC_CPUBackend<float>;
template class CMMC_CPUBackend<double>;


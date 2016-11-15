/*
 *  Created on: 26/nov/2010
 *      Author: flavio
 */

#include <set>
#include <sstream>

#include "VMMC_CPUBackend.h"
#include "../Utilities/Utils.h"
#include "../Interactions/DNAInteraction.h"
#include "../Particles/DNANucleotide.h"
#include "../Interactions/RNAInteraction2.h"
#include "../Interactions/DNA2Interaction.h"

template<typename number> VMMC_CPUBackend<number>::VMMC_CPUBackend() : MC_CPUBackend<number>() {
	//_op = NULL;
	_have_us = false;
	_netemps = 0;
	_etemps = NULL;
	_maxclust = 0;
	_reject_prelinks = false;
	_preserve_topology = false;
	_small_system = false;
	_last_move = MC_MOVE_TRANSLATION;
	_neighcells = NULL;
	_cells = NULL;
	_vmmc_heads = NULL;
	_U_ext = (number) 0.f;
	eijm = NULL;
	eijm_old = NULL;
	hbijm = NULL;
	hbijm_old = NULL;
	_default_weight = (number) 1.f;
	_safe_weights = true;
	_skip_hist_zeros = false;
	new_en3s = NULL;
	new_en5s = NULL;
	new_stn3s = NULL;
	new_stn5s = NULL;
	_vmmc_N_cells = 0;
	_equilibration_steps = 0;
}

template<typename number>
VMMC_CPUBackend<number>::~VMMC_CPUBackend() {
	// this is because otherwise the pointer to the force object gets freed
	// twice... maybe not the best way to do this, but oh well...
	if(this->_particles_old != NULL) {
		for (int i = 0; i < this->_N; i++) this->_particles_old[i]->N_ext_forces = 0;
	}
	_delete_cells();

	delete[] new_en3s;
	delete[] new_en5s;
	delete[] new_stn3s;
	delete[] new_stn5s;

	if (_netemps > 0)
		delete[] _etemps;

	if (_small_system) {
		if (eijm != NULL) {
			for (int k = 0; k < this->_N; k ++) delete[] eijm[k];
			delete[] eijm;
		}
		if (eijm_old != NULL) {
			for (int k = 0; k < this->_N; k ++) delete[] eijm_old[k];
			delete[] eijm_old;
		}
		if (hbijm != NULL) {
			for (int k = 0; k < this->_N; k ++) delete[] hbijm[k];
			delete[] hbijm;
		}
		if (hbijm_old != NULL) {
			for (int k = 0; k < this->_N; k ++) delete[] hbijm_old[k];
			delete[] hbijm_old;
		}
	}
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::init() {
	MC_CPUBackend<number>::init();

	// fix maxclust if evidently wrong
	if (_maxclust < 1) {
		OX_LOG(Logger::LOG_WARNING, "maxclust < 0, setting it to N = %i", this->_N);
		_maxclust = this->_N;
	}
	if (_maxclust > this->_N) {
		OX_LOG(Logger::LOG_WARNING, "maxclust > N does not make sense, setting it to N = %i", this->_N);
		_maxclust = this->_N;
	}

	if (_have_us) {
		_op.init_from_file(_op_file, this->_particles, this->_N);
		_w.init((const char *) _weights_file, &_op, _safe_weights, _default_weight);
		if (_reload_hist)
			_h.init(_init_hist_file, &_op, _etemps, _netemps);
		else
			_h.init(&_op, _etemps, _netemps);
		_h.set_simtemp(this->_T);
	}

	if (this->_delta[MC_MOVE_TRANSLATION] * sqrt(3) > this->_verlet_skin)
		throw oxDNAException("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");

	// setting the maximum displacement
	if (_preserve_topology) {
		_max_move_size = (number) 0.5;
		_max_move_size_sqr = _max_move_size * _max_move_size;
		OX_LOG(Logger::LOG_INFO, "Preserving topology; max_move_size = %lf...", _max_move_size);
	}
	else {
		_max_move_size = this->_box_side / 2. - 2. * this->_rcut - 0.2;
		_max_move_size_sqr = _max_move_size * _max_move_size;
		OX_LOG(Logger::LOG_INFO, "Not attempting to preserve topology; max_move_size = %g", _max_move_size);
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
	BaseParticle<number> * p, *q;
	for (int k = 0; k < this->_N; k ++) {
		p = this->_particles[k];
		if (p->n3 != P_VIRTUAL) {
			q = p->n3;
			epq = _particle_particle_bonded_interaction_n3_VMMC (p, q, &tmpf);
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

		for (int k = 0; k < this->_N; k ++) {
			for (int l = 0; l < k; l ++) {
				p = this->_particles[k];
				q = this->_particles[l];
				if (p->n3 != q && p->n5 != q) {
					eijm[k][l] = eijm[l][k] = eijm_old[k][l] = eijm_old[l][k] = _particle_particle_nonbonded_interaction_VMMC(p, q, &tmpf);
					hbijm[k][l] = hbijm[l][k] = hbijm_old[k][l] = hbijm_old[l][k] = (tmpf < HB_CUTOFF);
				}
			}
		}
	}

	_init_cells();

	this->_compute_energy();
	
	check_overlaps();

	if (this->_have_us) {
		this->_update_ops();
		check_ops();
	}
	
	if (this->_overlap == true)
		throw oxDNAException("There is an overlap in the initial configuration. Dying badly");
}


template<typename number>
void VMMC_CPUBackend<number>::get_settings(input_file & inp) {
	MC_CPUBackend<number>::get_settings(inp);
	int is_us, tmpi;

	std::string inter ("");
	if(getInputString(&inp, "interaction_type", inter, 0) == KEY_FOUND) {
		if(strncmp(inter.c_str(),"DNA2",512) != 0 && strncmp(inter.c_str(),"DNA2_nomesh",512) != 0 && strncmp(inter.c_str(),"RNA2",512) != 0 &&  strncmp(inter.c_str(), "DNA", 512) != 0 && strncmp(inter.c_str(), "DNA_nomesh", 512) != 0 && strncmp(inter.c_str(), "RNA", 512) != 0) throw oxDNAException("VMMC can be used only with DNA, DNA_nomesh, DNA2, DNA2_nomesh and RNA interactions");
	}

	if (getInputInt(&inp, "maxclust", &tmpi, 0) == KEY_FOUND) {
		_maxclust = tmpi;
		OX_LOG(Logger::LOG_INFO, "Using maxclust = %i", _maxclust);
	}

	if (getInputBoolAsInt(&inp, "small_system", &tmpi, 0) == KEY_FOUND) {
		if (tmpi > 0) {
			_small_system = true;
			OX_LOG(Logger::LOG_INFO, "Using algorithm N^2 for small system");
		}
		else {
			OX_LOG(Logger::LOG_INFO, "Using standard algorithm");
		}
	}

	if (getInputBoolAsInt(&inp, "preserve_topology", &tmpi, 0) == KEY_FOUND) {
		if (tmpi > 0) {
			_preserve_topology = true;
		}
		else {
			_preserve_topology = false;
		}
	}

	if (getInputBoolAsInt(&inp, "umbrella_sampling", &is_us, 0) != KEY_NOT_FOUND) {
		if (is_us > 0) {
			_have_us = true;
			getInputString(&inp, "op_file", _op_file, 1);
			getInputString(&inp, "weights_file", _weights_file, 1);
                       if (getInputString(&inp, "last_hist_file", _last_hist_file, 0) == KEY_NOT_FOUND) {
                               sprintf(_last_hist_file, "last_hist.dat"); OX_LOG(Logger::LOG_INFO, "Using default hist file %s", _last_hist_file);
                       }
                       if (getInputString(&inp, "traj_hist_file", _traj_hist_file, 0) == KEY_NOT_FOUND) {
                               sprintf(_traj_hist_file, "traj_hist.dat");
                               OX_LOG(Logger::LOG_INFO, "Using default traj hist file %s", _traj_hist_file);
                       }

                       // should we reload histograms?
                       if (getInputString(&inp, "init_hist_file", _init_hist_file, 0) == KEY_FOUND) {
                               OX_LOG(Logger::LOG_INFO, "Reloading histogram from %s", _init_hist_file);
                               this->_reload_hist = true;
                       } else {
                               this->_reload_hist = false;
                               FILE * tmpfile = fopen(_traj_hist_file, "w");
                               fclose(tmpfile);
                       }
			
			// wether to use unsafe weights
			if (getInputBoolAsInt (&inp, "safe_weights", &tmpi, 0) == KEY_FOUND) {
				_safe_weights = (tmpi > 0);
				if (! _safe_weights) {
					double tmpf;
					getInputDouble(&inp, "default_weight", &tmpf, 1);
					_default_weight = (number) tmpf;
					OX_LOG (Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) Using unsafe weights with default weight = %g", _default_weight);
				}
			}
			// wether to print out zero entries in traj_hist and last_hist 
			if (getInputBoolAsInt (&inp, "skip_hist_zeros", &tmpi, 0) == KEY_FOUND) {
				_skip_hist_zeros = tmpi > 0;
				if (_skip_hist_zeros) OX_LOG (Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) Skipping zero entries in traj_hist and last_hist files");
			} 
			
                       // should we extrapolate the histogram at different
                       // temperatures?
                       char tstring[512];
                       if (getInputString(&inp, "extrapolate_hist", tstring, 0) == KEY_FOUND) {
                               OX_LOG(Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) Extrapolating to temperatures .... %s", tstring);
                               char * aux, deg;
                               int c = 0, check;
                               double * tmpt;
                               tmpt = new double[100];
                               aux = strtok(tstring, ",");
                               while (aux != NULL) {
                                       //printf ("parsing %s\n", aux);
                                       check = sscanf(aux, "%lf %c", &(tmpt[c]), &deg);
                                       if (check < 1) {
                                               throw oxDNAException("(VMMC_CPUBackend.cpp) Unrecognizable line in extrapolate_hist");
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
                                                       throw oxDNAException("Unrecognizable temperature '%s' in extrapolate_hist", tmpt[c]);
                                                       break;
                                               }
                                       }
                                       c++;
                                       aux = strtok(NULL, ",");
                               }
                               if (c == 0) {
                                       throw oxDNAException("Nothing found in extrapolate_hist");
                               } else {
                                       char * tmps = (char *) malloc (100 * c * sizeof(char));
                                       char tmpss[512];
                                       sprintf (tmps, "%g ", tmpt[0]);
                                       for (int i = 1; i < c; i++) {
                                           sprintf (tmpss, "%g ", tmpt[i]); 
                                           strcat (tmps, tmpss);
                                       }
                                       OX_LOG (Logger::LOG_INFO, "Extrapolating to temperatures [in S. U.] %s", tmps);
                                       free (tmps);
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

	if (_have_us) {
		_h.read_interaction(inp);
	}

	if (getInputLLInt(&inp, "equilibration_steps", &_equilibration_steps, 0) == KEY_FOUND) {
		if (_equilibration_steps < 0) { 
			OX_LOG (Logger::LOG_WARNING, "found equilibration_steps < 0; setting equilibration steps = 0");
			_equilibration_steps = 0;
		}
		OX_LOG(Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) not recording histogram for steps < %lld", (long long unsigned)_equilibration_steps);
	}
}

// this function is just a wrapper that inverts p and q
template<typename number>
inline number VMMC_CPUBackend<number>::_particle_particle_bonded_interaction_n5_VMMC(BaseParticle<number> *p, BaseParticle<number> *q, number *stacking_en) {
	throw oxDNAException ("ERROR: called a function that should not be called; file %s, line %d", __FILE__, __LINE__);
	//return _particle_particle_bonded_interaction_n3_VMMC (q, p, stacking_en);
	//if (stacking_en != NULL) *stacking_en = 1.;
	//return this->_interaction->pair_interaction_bonded(q, p, NULL, false);
}

template<typename number>
inline number VMMC_CPUBackend<number>::_particle_particle_bonded_interaction_n3_VMMC(BaseParticle<number> *p, BaseParticle<number> *q, number *stacking_en) {
	BaseParticle<number> * tmp1, *tmp2, *tmp3, *tmp4;
	tmp1 = p->n3; tmp2 = p->n5; tmp3 = q->n3; tmp4 = q->n5;
	p->n3 = q;
	p->n5 = P_VIRTUAL;
	q->n3 = P_VIRTUAL;
	q->n5 = p;

	if (stacking_en != 0) *stacking_en = 0;

	LR_vector<number> r = q->pos - p->pos;

	// check for overlaps;
	LR_vector<number> rback = r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackr0;
        if(dynamic_cast<DNA2Interaction<number> *>(this->_interaction) != NULL) {
		rbackr0 = rback.module() - FENE_R0_OXDNA2;
	}
	else {
		rbackr0 = rback.module() - FENE_R0_OXDNA;
	}
	if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
		this->_overlap = true;
		p->n3 = tmp1;
		p->n5 = tmp2;
		q->n3 = tmp3;
		q->n5 = tmp4;
		return (number) (1.e6 * this->_T);
	}
	number energy = this->_interaction->pair_interaction_term (DNAInteraction<number>::BACKBONE, p, q, &r, false);

	energy += this->_interaction->pair_interaction_term (DNAInteraction<number>::BONDED_EXCLUDED_VOLUME, p, q, &r, false);

	number tmp_en = this->_interaction->pair_interaction_term (DNAInteraction<number>::STACKING, p, q, &r, false);
	energy += tmp_en;

	if (stacking_en != 0) *stacking_en = tmp_en;

	//if (true) {
		//number bef = this->_interaction->pair_interaction_bonded (p, q, NULL, false);
		//if (fabs ((bef - energy)) > 1.e-7 && bef < 10.)
			//printf ("%g %g  @@@\n ", bef, energy);
	//}

	p->n3 = tmp1;
	p->n5 = tmp2;
	q->n3 = tmp3;
	q->n5 = tmp4;

	return energy;
}

template<typename number>
inline number VMMC_CPUBackend<number>::_particle_particle_nonbonded_interaction_VMMC(BaseParticle<number> *p, BaseParticle<number> *q, number *H_energy) {

	//if (p->n3 == q || p->n5 == q || q->n3 == p || q->n5 == p)
	//	throw oxDNAException ("Called the function _particle_particel_nonbonded_interaction_VMMC with %s %d; %d, %d",__FILE__, __LINE__, p->index, q->index);

	if(H_energy != 0) *H_energy = (number) 0;

	LR_vector<number> r = q->pos.minimum_image (p->pos, this->_box_side);

	// early ejection
	if (r.norm() > this->_sqr_rcut) return (number) 0.f;

	number energy = this->_interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q, &r, false);

	if (H_energy != 0) *H_energy = energy;

	energy += this->_interaction->pair_interaction_term(DNAInteraction<number>::NONBONDED_EXCLUDED_VOLUME, p, q, &r, false);
	energy += this->_interaction->pair_interaction_term(DNAInteraction<number>::CROSS_STACKING, p, q, &r, false);

	// all interactions except DNA2Interaction use the DNAInteraction coaxial stacking
        if(dynamic_cast<DNA2Interaction<number> *>(this->_interaction) == NULL) energy += this->_interaction->pair_interaction_term(DNAInteraction<number>::COAXIAL_STACKING, p, q, &r, false);
	
	if(dynamic_cast<DNA2Interaction<number> *>(this->_interaction) != NULL)
	{
		energy += this->_interaction->pair_interaction_term(DNA2Interaction<number>::COAXIAL_STACKING, p, q, &r, false);
		energy += this->_interaction->pair_interaction_term(DNA2Interaction<number>::DEBYE_HUCKEL, p, q, &r, false);
	}
	else if(dynamic_cast<RNA2Interaction<number> *>(this->_interaction) != NULL)
	{
		
	        energy += this->_interaction->pair_interaction_term(RNA2Interaction<number>::DEBYE_HUCKEL, p, q, &r, false);
	}


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

template<typename number>
inline void VMMC_CPUBackend<number>::store_particle(BaseParticle<number> * src) {
	BaseParticle<number> *dst = this->_particles_old[src->index];

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->set_positions();
	dst->ext_potential = src->ext_potential;

	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::restore_particle (BaseParticle<number> * dst) {
	BaseParticle<number> *src = this->_particles_old[dst->index];

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->set_positions();
	dst->ext_potential = src->ext_potential;

	return;
}

template<typename number>
inline number VMMC_CPUBackend<number>::build_cluster_small (movestr<number> * moveptr, int maxsize, int * clust, int * size) {

	//printf ("before cycle...\n");
	//if (_have_us) check_ops();
	//printf ("passed\n");
	//printf ("\n\nfirst random: %g; %d @@@\n", drand48(), moveptr->seed);

	int nclust = 1;
	clust[0] = moveptr->seed;
	BaseParticle<number> * pp, *qq;
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

	//_check_metainfo();

	// CLUSTER GENERATION
	int k = 0;
	int neigh;
	pp = this->_particles[clust[0]];
	pp->inclust = true;

	//ppold = &(this->_particles_old[pp->index]);
	store_particle (pp);
	_move_particle(moveptr, pp, pp);

	while (k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = this->_particles[clust[k]];

		//printf ("recruiting from %i... @@@\n", pp->index);

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if (_preserve_topology) {
			if (pp->pos.sqr_distance (this->_particles_old[pp->index]->pos) > this->_max_move_size_sqr) {
				this->_dU = 0.;
				this->_dU_stack = 0.;
				* size = nclust;
				this->_overlap = true;
				return 0.;
			}
		}

		// trying to recruit bonded neighbors of pp
		if (pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
			if (!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC (pp, qq, &stack_temp);
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle (moveptr, qq, pp);

					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC (this->_particles_old[pp->index], qq);
					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes in the cluster
						//printf ("recruited %d ... @@@\n", qq->index);
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust ++;
					}
					else {
						//_r_move_particle(moveptr, qq);
						restore_particle (qq);
						prelinked_particles.insert(qq->index);
						//printf ("PRELINKED, NOT recruited %d - %g %g %g @@@\n", qq->index, E_old, E_pp_moved, E_qq_moved);
					}
				}
				else {
					new_en3s[pp->index] = new_en5s[qq->index] = E_pp_moved;
					new_stn3s[pp->index] = new_stn5s[qq->index] = stack_temp;
					//printf ("NOT PRELINKED, NOT recruited %d - %g %g @@@\n", qq->index, E_old, E_pp_moved);
				}
			}
		}

		// trying to recruit 5' neighbour of pp
		if (pp->n5 != P_VIRTUAL) {
			qq = pp->n5;
			if (!qq->inclust) {

				E_old = pp->en5;
				//E_pp_moved = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &stack_temp);
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC (qq, pp, &stack_temp);

				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle(moveptr, qq, pp);

					//E_qq_moved = _particle_particle_bonded_interaction_n5_VMMC (this->_particles_old[pp->index], qq);
					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC (qq, this->_particles_old[pp->index]);

					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes to cluster
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
						//printf ("recruited %d ... %g %g %g @@@\n", qq->index, E_old, E_pp_moved, E_qq_moved);
					}
					else {
						prelinked_particles.insert(qq->index);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
						//printf ("PRELINKED, NOT recruited %d ... @@@\n", qq->index);
					}
				}
				else {
					new_en5s[pp->index] = new_en3s[qq->index] = E_pp_moved;
					new_stn5s[pp->index] = new_stn3s[qq->index] = stack_temp;
					//printf ("NOT PRELINKED, NOT recruited %d ... @@@\n", qq->index);
					;
				}
			}
		}

		number tmpf = (number)0.;
		for (neigh = 0; neigh < this->_N; neigh ++) {
			qq = this->_particles[neigh]; //qq is my neighbor
			if ((!qq->inclust) && (pp->n3 != qq) && (pp->n5 != qq) && (qq != pp)) {

				E_old = eijm_old[pp->index][qq->index];

				if (E_old == (number)0.) continue;

				E_pp_moved = _particle_particle_nonbonded_interaction_VMMC (pp, qq, &tmpf);
				test1 = VMMC_link (E_pp_moved, E_old);
				if (test1 >  this->_next_rand ()) {
					store_particle (qq);
					_move_particle (moveptr, qq, pp);

					E_qq_moved = _particle_particle_nonbonded_interaction_VMMC (this->_particles_old[pp->index], qq);

					test2 = VMMC_link (E_qq_moved, E_old);
					if ((test2 / test1) > this->_next_rand()) {
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
						//printf ("recruited %d ... @@@\n", qq->index);
					}
					else {
						// prelinked;
						prelinked_particles.insert(qq->index);
						//_r_move_particle (moveptr, qq);
						restore_particle (qq);
						//printf ("PRELINKED, NOT recruited %d ... @@@\n", qq->index);
					}
				}
				else {
					//delta_E += E_pp_moved - E_old;
					eijm[pp->index][qq->index] = eijm[qq->index][pp->index] = E_pp_moved;
					hbijm[pp->index][qq->index] = hbijm[qq->index][pp->index] = (tmpf < HB_CUTOFF);
					//printf ("NOT PRELINKED, NOT recruited %d ... @@@\n", qq->index);
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
		//printf ("rejecting because of prelinks @@@");
		return (number) 0.;
	}

	// fix cells...
	for (int i = 0; i < nclust; i++) {
		int old_index, new_index;
		pp = this->_particles[clust[i]];
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
		pp = this->_particles[clust[i]];
		if(pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
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
			qq = pp->n5;
			if(qq->inclust == false) {
				//epq_new = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &tmpf_new);
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
			qq = this->_particles[neigh]; //qq is my neighbor
			if ((!qq->inclust) && (pp->n3 != qq) && (pp->n5 != qq) && (qq != pp)) {
				epq_old = eijm_old[pp->index][qq->index];

				if (epq_old == (number)0.) {
					// in this case we have not considered it in the
					// recruiting stage...
					epq_new = _particle_particle_nonbonded_interaction_VMMC (pp, qq, &tmpf_new);
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
	// N^2 algorithm, no problem, we don't need to implement anything

	pprime *= exp((1. / this->_T) * E_anomaly);

	//printf ("Returning %g @@@\n", pprime);

	this->_dU = delta_E;
	this->_dU_stack = delta_Est;

	return pprime;
}

// this function is the heart of the VMMC algorithm;this version uses
// cells to avoid the computation of O(N) interactions.
template<typename number>
inline number VMMC_CPUBackend<number>::build_cluster_cells (movestr<number> * moveptr, int maxsize, int * clust, int * size) {
	int nclust = 1;
	clust[0] = moveptr->seed;
	BaseParticle<number> * pp, *qq;
	number test1, test2;

	number pprime = 1.;

	number delta_E = 0;
	number delta_Est = 0;

	//check_ops();

	//this->_compute_energy();
	//this->_check_metainfo();
	//this->_check_old_metainfo();

	_reject_prelinks = false;

	set<int> prelinked_particles; //number of prelinked particles
	//set<base_pair, classcomp> poss_anomalies; //number of prelinked particles
	//set<base_pair, classcomp> poss_breaks; //number of prelinked particles
	set<base_pair, classcomp> prev_inter; // previous interactions

	number E_anomaly = 0;

	number E_qq_moved;
	number E_pp_moved;
	number E_old;
	number stack_temp;
	number H_temp;

	// CLUSTER GENERATION
	int k = 0;
	int icell, neigh;
	pp = this->_particles[clust[0]];
	pp->inclust = true;

	store_particle (pp);
	_move_particle(moveptr, pp, pp);

	assert (this->_overlap == false);

	// how far away am I recruiting?
	LR_vector<number> how_far (0., 0., 0.);

	while (k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = this->_particles[clust[k]];

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if (_preserve_topology) {
			if (pp->pos.sqr_distance (this->_particles_old[pp->index]->pos) > this->_max_move_size_sqr) {
				this->_dU = 0.;
				this->_dU_stack = 0.;
				* size = nclust;
				this->_overlap = true;
				return 0.;
			}
		}

		// trying to recruit bonded neighbors of pp
		if (pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
			if (!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC (pp, qq, &stack_temp);
				test1 = VMMC_link (E_pp_moved, E_old);
				if (this->_overlap || test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle(moveptr, qq, pp);

					// in case E_pp_moved created an overlap
					this->_overlap = false;

					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC (this->_particles_old[pp->index], qq);
					//_move_particle(moveptr, pp);

					test2 = VMMC_link (E_qq_moved, E_old);
					if (this->_overlap || (test2 / test1) > this->_next_rand()) {
						//we did full_link, qq goes in the cluster

						// in case E_qq_moved created an overlap
						this->_overlap = false;

						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust ++;
					}
					else {
						assert (this->_overlap == false);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
						prelinked_particles.insert(qq->index);
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
			qq = pp->n5;
			if (!qq->inclust) {

				E_old = pp->en5;
				//E_pp_moved = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &stack_temp);
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC (qq, pp, &stack_temp);

				test1 = VMMC_link (E_pp_moved, E_old);
				if (this->_overlap || test1 > this->_next_rand()) {
					// prelink successful
					store_particle (qq);
					_move_particle(moveptr, qq, pp);

					// in case we have recruited because of an overlap
					this->_overlap = false;

					//_r_move_particle(moveptr, pp);
					//E_qq_moved = _particle_particle_bonded_interaction_n5_VMMC (pp, qq);
					//_move_particle(moveptr, pp);
					//E_qq_moved = _particle_particle_bonded_interaction_n5_VMMC (this->_particles_old[pp->index], qq);
					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC (qq, this->_particles_old[pp->index]);

					test2 = VMMC_link (E_qq_moved, E_old);
					if (this->_overlap || (test2 / test1) > this->_next_rand()) {
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
			//icell = cell_neighbours (_cells[pp->index], c);
			//assert (icell == _neighcells[_cells[pp->index]][c]);
			neigh = _vmmc_heads[icell];
			while (neigh != P_INVALID) {
				qq = this->_particles[neigh]; //qq is my neighbor

				if (pp->n3 == qq || pp->n5 == qq) {
					neigh = qq->next_particle;
					continue;
				}

				if (qq->inclust == false) {
					E_old = _particle_particle_nonbonded_interaction_VMMC (this->_particles_old[pp->index], qq, &H_temp);

					if (E_old == (number)0.) {
						neigh = qq->next_particle;
						continue;
					}

					E_pp_moved = _particle_particle_nonbonded_interaction_VMMC (pp, qq);

					test1 = VMMC_link (E_pp_moved, E_old);
					if (test1 >  this->_next_rand ()) {
						store_particle (qq);
						_move_particle (moveptr, qq, pp);

						//_r_move_particle (moveptr, pp);
						//E_qq_moved = _particle_particle_nonbonded_interaction_VMMC (pp, qq);
						//_move_particle (moveptr, pp);
						E_qq_moved = _particle_particle_nonbonded_interaction_VMMC (this->_particles_old[pp->index], qq);

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
						if (fabs (E_old) > 0.) {
							// we store the possible interaction to account for later
							store_particle (qq);
							prev_inter.insert ((pp->index > qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						}
					}
				}
				neigh = qq->next_particle;
			}
		}
		k ++;
	}
	*size = nclust;

	if (this->_overlap) printf ("cera anche qui\n");

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
		pp = this->_particles[clust[i]];
		old_index = _cells[pp->index];
		new_index = _get_cell_index(pp->pos);
		if (new_index != old_index) {
			_fix_list (pp->index, old_index, new_index);
		}
	}

	delta_E = delta_Est = E_anomaly = 0.;
	number tmpf_old, tmpf_new, epq_new, epq_old;
	bool h_new, h_old;
	for (int i = 0; i < nclust; i++) {
		pp = this->_particles[clust[i]];
		if(pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
			if(qq->inclust == false) {
				epq_new = _particle_particle_bonded_interaction_n3_VMMC (pp, qq, &tmpf_new);

				assert (this->_overlap == false);
				delta_E += epq_new - pp->en3;
				delta_Est += tmpf_new - pp->esn3;

				new_en3s[pp->index] = new_en5s[qq->index] = epq_new;
				new_stn3s[pp->index] = new_stn5s[qq->index] = tmpf_new;
			}
		}

		if(pp->n5 != P_VIRTUAL) {
			qq = pp->n5;
			if(qq->inclust == false) {
				if (this->_overlap) printf ("cera prima\n");
				epq_new = _particle_particle_bonded_interaction_n3_VMMC (qq, pp, &tmpf_new);
				//epq_new = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &tmpf_new);

				assert (this->_overlap == false);

				delta_E += epq_new - pp->en5;
				delta_Est += tmpf_new - pp->esn5;

				new_en5s[pp->index] = new_en3s[qq->index] = epq_new;
				new_stn5s[pp->index] = new_stn3s[qq->index] = tmpf_new;
			}
		}

		for (int c = 0; c < 27; c ++) {
			icell = _neighcells[_cells[pp->index]][c];
			//icell = cell_neighbours (_cells[pp->index], c);
			//assert (icell == _neighcells[_cells[pp->index]][c]);
			neigh = _vmmc_heads[icell];
			while (neigh != P_INVALID) {
				qq = this->_particles[neigh];

				if (pp->n3 == qq || pp->n5 == qq) {
					neigh = qq->next_particle;
					continue;
				}

				if (qq->inclust == false) {
					//_r_move_particle (moveptr, pp);
					//epq_old = _particle_particle_nonbonded_interaction_VMMC (pp, qq, &tmpf_old);
					//_move_particle (moveptr, pp);
					epq_old = _particle_particle_nonbonded_interaction_VMMC (this->_particles_old[pp->index], qq, &tmpf_old);
					epq_new = _particle_particle_nonbonded_interaction_VMMC (pp, qq, &tmpf_new);

					delta_E += epq_new - epq_old;

					// we have considered this interaction, so we remove it from the list
					if (fabs (epq_old) > 0.) prev_inter.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));

					// check for anomaly of second kind;
					if (epq_old == 0. && epq_new > 0.) {
						// we have just created an overlap where there
						// was no interaction
						E_anomaly -= epq_new;
					}

					// check for anomaly of first kind
					if (epq_old > 0.) {
						if (epq_new == 0.) {
							E_anomaly += epq_old;
						}
					}
					// fix h_bonding...
					if (_have_us) {
						h_new = tmpf_new < HB_CUTOFF;
						h_old = tmpf_old < HB_CUTOFF;
						//poss_breaks.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
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
				neigh = qq->next_particle;
			}
		}
	}

	// now we treat the case where we moved A LOT; in this case,
	// a particle pair that had an interaction different from 0. before the move
	// can now be such a high distance in between that that pair is not considered
	// in the computation of delta_E above. In this case, we use the information
	// about the previously interacting pairs to account for those contributions to
	// delta_E. Also, this allows us to update the anomalies of the fist kind (removal
	// of positive interactions) and also to update the order parameter, in case we
	// have an hydrogen bond that was present before the move in between nucleotides
	// that are now far apart
	set<base_pair>::iterator it;
	for (it = prev_inter.begin(); it != prev_inter.end(); it ++) {
		number tmpf;
		pp = this->_particles[(*it).first];
		qq = this->_particles[(*it).second];
		if (!(pp->inclust && qq->inclust)) {
			epq_old = _particle_particle_nonbonded_interaction_VMMC (this->_particles_old[pp->index], this->_particles_old[qq->index], &tmpf);
			delta_E -= epq_old;
			if (epq_old > 0.) E_anomaly += epq_old;
			// if we get to this stage, epq_new has to be 0, so the hydrogen bond
			// is for sure not present anymore
			if (_have_us) {
				if (tmpf < HB_CUTOFF) {
					_op.remove_hb (pp->index, qq->index);
				}
			}
		}
	}

	pprime *= exp((1. / this->_T) * E_anomaly);

	this->_dU = delta_E;
	this->_dU_stack = delta_Est;

	return pprime;
}

/* // old version
template<typename number>
inline void VMMC_CPUBackend<number>::_move_particle(movestr<number> * moveptr, BaseParticle< number> *q) {
	if (moveptr->type == MC_MOVE_TRANSLATION) {
		q->pos += moveptr->t;
	}
	else if (moveptr->type == MC_MOVE_ROTATION) {
		//in this case, the translation vector is the point around which we
		//rotate
		LR_vector<number> dr, drp;
		if (this->_particles[moveptr->seed]->strand_id == q->strand_id) {
			dr = q->pos - moveptr->t;
		//if (dr.norm() > 29. * 29.) {
		//		printf("cavolo1... %g %g %g %d %d\n", dr.x, dr.y, dr.z, moveptr->seed, q->index);
		//	}
		}
		else {
			dr = q->pos.minimum_image(moveptr->t, this->_box_side);
		//	if (dr.norm() > 29. * 29.) {
		//		printf("cavolo2... %g %g %g %d %d\n", dr.x, dr.y, dr.z, moveptr->seed, q->index);
		//	}
		}

		if (q->index == 122 || q->index == 121) printf ("rotating %d... ", q->index);
		if (q->index == 122 || q->index == 121) printf ("pos  = np.array ([%g, %g, %g]\n", q->pos.x, q->pos.y, q->pos.z);
		if (q->index == 122 || q->index == 121) printf ("da:  = np.array ([%g, %g, %g]\n", moveptr->t.x, moveptr->t.y, moveptr->t.z);
		if (q->index == 122 || q->index == 121) printf ("dr:  = np.array ([%g, %g, %g]\n", dr.x, dr.y, dr.z); 
		drp = moveptr->R * dr;
		if (q->index == 122 || q->index == 121) printf ("drp: = np.array ([%g, %g, %g]\n", drp.x, drp.y, drp.z); 
		q->pos += (drp - dr); // accounting for PBC
		q->orientation = moveptr->R * q->orientation;
		//q->orientation.orthonormalize(); // NORMALIZZATo
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();
		if (q->index == 122 || q->index == 121) printf ("fin: = np.array ([%g, %g, %g]\n", q->pos.x, q->pos.y, q->pos.z);
	}
	else {
		;
	}

	return;
}
*/

template<typename number>
inline void VMMC_CPUBackend<number>::_move_particle(movestr<number> * moveptr, BaseParticle< number> *q, BaseParticle<number> *p) {
	if (moveptr->type == MC_MOVE_TRANSLATION) {
		q->pos += moveptr->t;
	}
	else if (moveptr->type == MC_MOVE_ROTATION) {
		//in this case, the translation vector is the point around which we rotate
		BaseParticle<number> * p_old = this->_particles_old[p->index];
		LR_vector<number> dr, drp;
		if (p->strand_id == q->strand_id) { dr = q->pos - p_old->pos; }
		else { dr = q->pos.minimum_image(p_old->pos, this->_box_side); }

		// we move around the backbone site
		//dr += moveptr->t;

		drp = moveptr->R * dr;
		//q->pos += (drp - dr); // accounting for PBC
		q->pos = p->pos + drp;
		q->orientation = moveptr->R * q->orientation;
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();
	}
	else {
		;
	}

	return;
}

template<typename number>
inline void VMMC_CPUBackend<number>::_fix_list (int p_index, int oldcell, int newcell) {
	int j, jold;

	// remove p_index from its old cell
	//printf ("## %i %i %i\n", oldcell, _vmmc_N_cells, _cells[p_index]);
	j = _vmmc_heads[oldcell];
	jold = P_INVALID;
	assert (j != P_INVALID);
	while (j != p_index) {
		jold = j;
		j = this->_particles[j]->next_particle;
	}

	//printf ("j, jold: %i, %i\n", j, jold);
	//printf("no\n");
	assert (j != jold);
	assert (j != P_INVALID);

	if (jold != P_INVALID) {
		this->_particles[jold]->next_particle = this->_particles[p_index]->next_particle;
	}
	else {
		_vmmc_heads[oldcell] = this->_particles[p_index]->next_particle;
	}

	// the old was full, since it contained the particle; is it now empty?
	if (_vmmc_heads[oldcell] == P_INVALID) {
		//printf ("now empty %i\n", oldcell);
		delete [] _neighcells[oldcell];
	}
	if (_vmmc_heads[newcell] == P_INVALID) {
		//printf ("now filled %i\n", newcell);
		_neighcells[newcell] = new int [27];
		int ind[3], loop_ind[3], nneigh = 0;
		ind[0] = newcell % _vmmc_N_cells_side;
		ind[1] = (newcell / _vmmc_N_cells_side) % _vmmc_N_cells_side;
		ind[2] = newcell / (_vmmc_N_cells_side * _vmmc_N_cells_side);
		for(int jj = -1; jj < 2; jj++) {
			loop_ind[0] = (ind[0] + jj + _vmmc_N_cells_side) % _vmmc_N_cells_side;
			for(int kk = -1; kk < 2; kk++) {
				loop_ind[1] = (ind[1] + kk + _vmmc_N_cells_side) % _vmmc_N_cells_side;
				for(int ll = -1; ll < 2; ll++) {
					loop_ind[2] = (ind[2] + ll + _vmmc_N_cells_side) % _vmmc_N_cells_side;
					int loop_index = (loop_ind[2] * _vmmc_N_cells_side + loop_ind[1]) * _vmmc_N_cells_side + loop_ind[0];
					_neighcells[newcell][nneigh] = loop_index;
					nneigh ++;
				}
			}
		}
	}

	// add it to the new cell
	this->_particles[p_index]->next_particle = _vmmc_heads[newcell];
	_vmmc_heads[newcell] = p_index;

	_cells[p_index] = newcell;

	return;
}


template<typename number>
void VMMC_CPUBackend<number>::sim_step(llint curr_step) {

	this->_mytimer->resume();
	this->_timer_move->resume();

	//srand48(curr_step);
	/*
	  {
		unsigned short seednow[3];
		Utils::get_seed((unsigned short *)seednow);
		printf ("%hu %hu %hu\n", seednow[0], seednow[1], seednow[2]);
	  }*/

	LR_vector<number> tmp;
	//BaseParticle<number> *pold;

	//printf ("checking metainfo at the beginning...\n");
	//_check_metainfo();
	//printf ("checking metainfo at the beginning: PASSED\n");

	//	_compute_energy();
	//	printf("%lf %lf\n", this->_U/this->_N, this->_U_hydr);
	//	exit(1);

	//get_time(&this->_timer, 0);

	int * clust, nclust;
	clust = new int[this->_N];

	double oldweight, weight;
	int windex, oldwindex;
	oldweight = weight = 1.;
	if (_have_us) oldweight = _w.get_weight(_op.get_all_states(), &oldwindex);

	/*
	// normalisation factor for 1/_maxclust
	number norm = log(_maxclust) + (number) 0.5772156649 + (number) (1. / 2.)
			/ (number) _maxclust + (number) (1. / 12.) / _maxclust / _maxclust;
	number logmaxclustplusone = (number) log((number) _maxclust + 1.);
	number M = (1. / norm) / (1. / (logmaxclustplusone * 2.001));
	*/
	
	// set the potential due to external forces
	_U_ext = (number) 0.f;
	for (int k = 0; k < this->_N; k ++) {
		BaseParticle<number> *p = this->_particles[k];
		p->set_ext_potential(curr_step, this->_box_side);
		_U_ext += p->ext_potential;
	}

	for (int i = 0; i < this->_N; i++) {
		if (_have_us) _op.store();
		this->_dU_stack = 0.;
		//printf ("\n##A %lf %lf \n", this->_U_stack, this->_dU_stack);

		// check of ext // works with TWO traps, not with one
		//number ov_c = (number) 0;
		//for (int l = 0; l < this->_N; l ++) {
		//	BaseParticle<number> * pp = &(this->_particles[l]);
		//	pp->set_ext_potential(curr_step);
		//	ov_c += pp->ext_potential;
		//}

		// seed particle;
		int pi = (int) (drand48() * this->_N);
		BaseParticle<number> *p = this->_particles[pi];

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
			//move.t = this->_particles[move.seed]->int_centers[DNANucleotide<number>::BACK] + this->_particles[move.seed]->pos;
			move.t = this->_particles[move.seed]->int_centers[DNANucleotide<number>::BACK];
			if (fabs((move.t * move.t)) > 0.5) printf("caca");
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
		if (this->_overlap == false && pprime > 0.) {
			for (int l = 0; l < nclust; l++) {
				p = this->_particles[clust[l]];
				delta_E_ext += - p->ext_potential;
				p->set_ext_potential(curr_step, this->_box_side);
				delta_E_ext += + p->ext_potential;
			}
			pprime *= exp(-(1. / this->_T) * delta_E_ext);
		}

		_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

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
			if (nclust <= _maxclust) this->_accepted[_last_move]++;
			//if (!_reject_prelinks) this->_accepted[0]++;
			this->_U += this->_dU;
			this->_U_stack += this->_dU_stack;
			_U_ext += delta_E_ext;

			oldweight = weight; // if (!_have_us) oldweight = weight = 1.;
			oldwindex = windex; // if (!_have_us) oldweight = weight = 1.;

			for (int l = 0; l < nclust; l ++) {
				BaseParticle<number> * pp, * qq;
				pp = this->_particles[clust[l]];
				if (pp->n3 != P_VIRTUAL) {
					qq = pp->n3;
					if (qq->inclust == false) {
						pp->en3 = qq->en5 = new_en3s[clust[l]];
						pp->esn3 = qq->esn5 = new_stn3s[clust[l]];
					}
				}
				if (pp->n5 != P_VIRTUAL) {
					qq = pp->n5;
					if (qq->inclust == false) {
						pp->en5 = qq->en3 = new_en5s[clust[l]];
						pp->esn5 = qq->esn3 = new_stn5s[clust[l]];
					}
				}

				if (_small_system) {
					for (int c = 0; c < this->_N; c ++) {
						qq = this->_particles[c];
						if (pp != qq && pp->n3 != qq && pp->n5 != qq && qq->inclust == false) {
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
				BaseParticle<number> * pp;
				pp = this->_particles[clust[l]];
				//_r_move_particle (&move, pp);
				restore_particle (pp);
				old_index = _cells[pp->index];
				new_index = _get_cell_index(pp->pos);
				if (new_index != old_index) {
					_fix_list (pp->index, old_index, new_index);
				}
				pp->set_ext_potential(curr_step, this->_box_side);
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
			BaseParticle<number> * pp;
			pp = &(this->_particles[l]);
			c_ext += pp->ext_potential;
			pp->set_ext_potential(curr_step);
			c_ext_fs += pp->ext_potential;
		}
		if (fabs (c_ext - c_ext_fs) > 1.e-6) {
			fprintf (stderr, "%g %g -- beh\n", c_ext, c_ext_fs);
		}
		// check ext potential done*/

		// add to the histogram
		if (_have_us && curr_step > _equilibration_steps) _h.add(oldwindex, oldweight, this->_U, this->_U_stack, _U_ext);

		// reset the inclust property to the particles
		for (int k = 0; k < nclust; k++) {
			this->_particles[clust[k]]->inclust = false;
		}
	}

	//check_ops();

	delete[] clust;
	//delete[] tainted;
	
	// check energy for percolation
	if (curr_step % (llint)this->_check_energy_every == 1) {
		//printf ("checking energy for percolation..\n");
		number U_from_tally = this->_U;
		_compute_energy();
		if ((this->_U - U_from_tally) > 1.e-4) throw oxDNAException ("(VMMC_CPUBackend) Accumulated Energy (%g) and Energy computed from scratch (%g) don't match. Possibly percolating clusters. Your box is too small", U_from_tally, this->_U);
		//printf ("all ok (%g %g)... \n", U_from_tally, this->_U);
	}

	//get_time(&this->_timer, 1);
	//process_times(&this->_timer);
	this->_timer_move->pause();
	this->_mytimer->pause();
}

template<typename number>
void VMMC_CPUBackend<number>::check_ops() {
	if (!_have_us) return;
	//printf ("checking OP...\n");
	assert (_have_us);

	int * state;

	//state = (int *) malloc(_op.get_hb_parameters_count() * sizeof(int));
	//memcpy(state, _op.get_hb_states(), _op.get_hb_parameters_count() * sizeof(int));
	state = (int *) malloc(_op.get_all_parameters_count() * sizeof(int));
	memcpy(state, _op.get_all_states(), _op.get_all_parameters_count() * sizeof(int));
	_op.reset();

	int i, j;
	BaseParticle<number> *p, *q;
	number hpq;
	for (i = 0; i < this->_N; i++) {
		p = this->_particles[i];
		for (j = 0; j < i; j ++) {
			q = this->_particles[j];
			if (p->n3 != q && p->n5 != q) {
				_particle_particle_nonbonded_interaction_VMMC (p, q, &hpq);
				if (hpq < HB_CUTOFF) {
					_op.add_hb (p->index, q->index);
				}
			}
		}
	}
	_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

	int * new_state = _op.get_all_states ();
	int check = 0;

	for (i = 0; i < _op.get_all_parameters_count(); i++) {
		if (state[i] != new_state[i])
			printf("%d should be %d \n", state[i], new_state[i]);
		check += abs(new_state[i] - state[i]);
	}

	/*
	int * new_state = _op.get_hb_states();
	int check = 0;
	for (i = 0; i < _op.get_hb_parameters_count(); i++) {
		if (state[i] != new_state[i])
			printf("%d should be %d \n", state[i], new_state[i]);
		check += abs(new_state[i] - state[i]);
	}*/

	if (check != 0) {
		printf ("CASINO\n");
		abort();
	}

	free(state);
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_update_ops() {

	assert (_have_us);

	_op.reset();

	// hydrogen bonding
	int i, j, c;
	BaseParticle<number> *p, *q;
	number hpq;
	for (i = 0; i < this->_N; i++) {
		p = this->_particles[i];
		for (c = 0; c < 27; c++ ) {
			j = _vmmc_heads[_neighcells[_cells[p->index]][c]];
			while (j != P_INVALID) {
				q = this->_particles[j];
				if (j < i && p->n3 != q && p->n5 != q) {
					_particle_particle_nonbonded_interaction_VMMC (p, q, &hpq);
					if (hpq < HB_CUTOFF) {
						_op.add_hb (i, j);
					}
				}
				j = q->next_particle;
			}
		}
	}

	// distances
	_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

	//exit(-1);
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::_print_pos(int id) {
	BaseParticle<number> * p;
	p = this->_particles[id];
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
inline void VMMC_CPUBackend<number>::check_overlaps() {
	int i, noverlaps;
	//number epq;
	BaseParticle<number> *p, *q;

	noverlaps = 0;
	for (i = 0; i < this->_N; i ++) {
		p = this->_particles[i];
		if (p->n3 != P_VIRTUAL) {
			this->_overlap = false;
			q = p->n3;
			_particle_particle_bonded_interaction_n3_VMMC (p, q);
			if (this->_overlap) {
				noverlaps ++;
				LR_vector<number> rbb, r;
				r = p->pos - q->pos;
				rbb = r + p->int_centers[DNANucleotide<number>::BACK] - q->int_centers[DNANucleotide<number>::BACK];
				printf ("### overlap %i and %i (%g, %g)\n", p->index, q->index, r.module(), rbb.module());
				this->_overlap = false;
			}
		}
	}
	assert (noverlaps == 0);
	if (noverlaps > 0) abort();
}

template<typename number>
char * VMMC_CPUBackend<number>::get_op_state_str() {
	if (_have_us) {
		int * state = _op.get_all_states();
		char * aux;
		aux = (char *) _state_str;
		for (int i = 0; i < _op.get_all_parameters_count(); i++) {
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
void VMMC_CPUBackend<number>::print_conf(llint curr_step, bool reduced,	bool only_last) {
	SimBackend<number>::print_conf(curr_step, reduced, only_last);
	if (_have_us) {
		if (!only_last) _h.print_to_file(_traj_hist_file, curr_step, false, _skip_hist_zeros);
		_h.print_to_file(_last_hist_file, curr_step, true, _skip_hist_zeros);
	}
}

template<typename number>
void VMMC_CPUBackend<number>::print_conf(llint curr_step, bool only_last) {
	SimBackend<number>::print_conf(curr_step, only_last);
	if (_have_us) {
		if (!only_last)
			this->_h.print_to_file(_traj_hist_file, curr_step, false, _skip_hist_zeros);
		_h.print_to_file(_last_hist_file, curr_step, true, _skip_hist_zeros);
	}
}

template<typename number>
void VMMC_CPUBackend<number>::_compute_energy () {
	// Since this function is called by MC_CPUBackend::init() but it uses cells initialized
	// by VMMC_CPUBackend::init(), we do nothing if it's called too early
	
	if(_vmmc_heads == NULL) return;
	BaseParticle<number> * p, * q;
	this->_overlap = false;
	number res = (number) 0;
	number dres, tmpf;

	this->_U = this->_U_hydr = (number) 0;

	for (int i = 0; i < this->_N; i ++) {
		p = this->_particles[i];
		if (p->n3 != P_VIRTUAL) {
			q = p->n3;
			dres = _particle_particle_bonded_interaction_n3_VMMC (p, q);
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
			while (j != P_INVALID) {
				q = this->_particles[j];
				if (p->n3 != q && p->n5 != q && p->index < q->index) {
					dres = _particle_particle_nonbonded_interaction_VMMC (p, q, &tmpf);
					this->_U += dres;
					this->_U_hydr += tmpf;
				}
				j = q->next_particle;
			}
		}
	}

	if (this->_overlap) {
	    throw oxDNAException ("overlap found. Aborting..\n");
	}

}

template<typename number>
void VMMC_CPUBackend<number>::_init_cells() {
	//_vmmc_N_cells_side = (int) floor(this->_box_side / sqrt(this->_sqr_rverlet));
	_vmmc_N_cells_side = (int) (floor(this->_box_side / this->_rcut) + 0.1);
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
	while(_vmmc_N_cells_side > ceil(750.) && _vmmc_N_cells_side > 3) {
		_vmmc_N_cells_side--;
	}

	if(_vmmc_N_cells_side < 3) throw oxDNAException("N_cells_side (%d) must be > 2. You're gonna need a bigger box.", _vmmc_N_cells_side);

	_vmmc_N_cells = _vmmc_N_cells_side * _vmmc_N_cells_side * _vmmc_N_cells_side;

	_vmmc_heads = new int[_vmmc_N_cells];
	_cells = new int[this->_N];
	_neighcells = new int * [_vmmc_N_cells];
	//for (int i = 0; i < _vmmc_N_cells; i ++)
	// _neighcells[i] = new int[27];

	/*
	for(int i = 0; i < _vmmc_N_cells; i++) {
		_cells[i][0] = i % _vmmc_N_cells_side;
		_cells[i][1] = i / _vmmc_N_cells_side;
		_cells[i][2] = i / (_vmmc_N_cells_side * _vmmc_N_cells_side);
		fprintf (stderr, "CELLE: %i -> (%i, %i, %i)\n", i, _cells[i][0], _cells[i][1], _cells[i][2]);
	}*/

	//fprintf (stderr, "VMMC CELL INFO: N_cells=%i, box_side=%g, N_cells_side=%i, r_cut = %g\n", _vmmc_N_cells, this->_box_side, _vmmc_N_cells_side, this->_rcut);
	/*
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
				}
			}
		}
		assert (nneigh == 27);
	}*/
	//fprintf (stderr, "cells initialised...\n");

	for(int i = 0; i < _vmmc_N_cells; i++)
	  _vmmc_heads[i] = P_INVALID;

	for(int i = 0; i < this->_N; i++) {
		this->_particles[i]->next_particle = P_INVALID;
	}

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		int cell_index = _get_cell_index(p->pos);
		int old_head = _vmmc_heads[cell_index];
		_vmmc_heads[cell_index] = i;
		_cells[i] = cell_index;
		assert (cell_index < _vmmc_N_cells);
		p->next_particle = old_head;
	}

	for (int i = 0; i < _vmmc_N_cells; i ++) {
		if (_vmmc_heads[i] != P_INVALID) {
			_neighcells[i] = new int [27];
			int ind[3], loop_ind[3], nneigh = 0;
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
					}
				}
			}
		}
	}

	// check for the cells
	for (int i = 0; i < _vmmc_N_cells; i ++) {
		//fprintf (stderr, "## i %i\n", i);
		int j = _vmmc_heads[i];
		//fprintf (stderr, "## j %i\n", j);
		if (j != P_INVALID) {
			// fprintf (stderr, "cell %5i: %i ", i, j);
			j = this->_particles[j]->next_particle;
			while (j != P_INVALID) {
				//fprintf (stderr, "%i ", j);
				j = this->_particles[j]->next_particle;
			}
			// fprintf (stderr, "\n");
		}
	}
	return;
}


template<typename number>
void VMMC_CPUBackend<number>::_delete_cells() {
	//for (int i = 0; i < this->_vmmc_N_cells; i ++) delete[] _neighcells[i];
	for (int i = 0; i < this->_vmmc_N_cells; i ++) {
		if (_vmmc_heads[i] != P_INVALID) {
			delete[] _neighcells[i];
		}
	}
	delete [] _vmmc_heads;
	delete [] _cells;
	delete [] _neighcells;
	return;
}

template<typename number>
void VMMC_CPUBackend<number>::fix_diffusion() {
	
	// fix diffusion can sometimes change the value of the order paramer by changing the
	// orientations/coordinates particles that were barely above/below a cutoff.
	// To avoid this, we store the system before fix_diffusion and restore it 
	// if the order parameter value has changed.

	for (int i = 0; i < this->_N; i ++) store_particle (this->_particles[i]);

	// check of order parameter
	int * state = (int *) malloc(_op.get_all_parameters_count() * sizeof(int));
	memcpy(state, _op.get_all_states(), _op.get_all_parameters_count() * sizeof(int));
	_op.reset();

	SimBackend<number>::fix_diffusion();

	_op.reset();
	for (int i = 0; i < this->_N; i++) {
		BaseParticle<number> * p = this->_particles[i];
		for (int j = 0; j < i; j ++) {
			BaseParticle<number> * q = this->_particles[j];
			if (p->n3 != q && p->n5 != q) {
				number hpq;
				_particle_particle_nonbonded_interaction_VMMC (p, q, &hpq);
				if (hpq < HB_CUTOFF) _op.add_hb (p->index, q->index);
			}
		}
	}
	_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

	int * new_state = _op.get_all_states ();

	int check = 0;
	for (int i = 0; i < _op.get_all_parameters_count(); i++) check += abs(new_state[i] - state[i]);

	if (check != 0) {
		OX_LOG (Logger::LOG_DEBUG, "(VMMC_CPUBackend) fix_diffusion() changed the value of the order parameter. Restoring simulation status before fix_diffusion()");
		for (int i = 0; i < this->_N; i ++) restore_particle (this->_particles[i]);
	}
}

template<typename number>
void VMMC_CPUBackend<number>::print_observables(llint curr_step) {
	this->_lists->global_update(true);
	this->_backend_info += get_op_state_str();
	MCBackend<number>::print_observables(curr_step);
	//if ((curr_step % (10 * this->_N)) == 0) this->fix_diffusion();
}

template<>
inline int VMMC_CPUBackend<float>::_get_cell_index(const LR_vector<float> &pos) {
	int res = (int) ((pos.x / _box_side - floor(pos.x / _box_side)) * (1.f - FLT_EPSILON) * _vmmc_N_cells_side);
	res += _vmmc_N_cells_side * ((int) ((pos.y / _box_side - floor(pos.y / _box_side)) * (1.f - FLT_EPSILON) * _vmmc_N_cells_side));
	res += _vmmc_N_cells_side * _vmmc_N_cells_side *
		((int) ((pos.z / _box_side - floor(pos.z / _box_side)) * (1.f - FLT_EPSILON) * _vmmc_N_cells_side));
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


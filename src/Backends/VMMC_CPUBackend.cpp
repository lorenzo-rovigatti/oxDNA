/*
 *  Created on: 26/nov/2010
 *      Author: flavio
 */

#include "VMMC_CPUBackend.h"
#include "../Utilities/Utils.h"
#include "../Interactions/RNAInteraction2.h"
#include "../Interactions/DNA2Interaction.h"
#include "../Interactions/DRHInteraction.h"

#include <set>
#include <sstream>
#include <limits>

VMMC_CPUBackend::VMMC_CPUBackend() :
				MC_CPUBackend() {
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
	_vmmc_N_cells = 0;
	_vmmc_box_side = -1.;
	_equilibration_steps = 0;
	_vmmc_N_cells_side = -1;
	_reload_hist = false;
	_just_updated_lists = false;

	_dU = 0.;
	_U_stack = 0.;
	_dU_stack = 0.;

	_max_move_size = 0.5;
	_max_move_size_sqr = SQR(_max_move_size);
}

VMMC_CPUBackend::~VMMC_CPUBackend() {
	_delete_cells();

	if(_netemps > 0) {
		delete[] _etemps;
	}

	if(_small_system) {
		if(eijm != NULL) {
			for(int k = 0; k < N(); k++)
				delete[] eijm[k];
			delete[] eijm;
		}
		if(eijm_old != NULL) {
			for(int k = 0; k < N(); k++)
				delete[] eijm_old[k];
			delete[] eijm_old;
		}
		if(hbijm != NULL) {
			for(int k = 0; k < N(); k++)
				delete[] hbijm[k];
			delete[] hbijm;
		}
		if(hbijm_old != NULL) {
			for(int k = 0; k < N(); k++)
				delete[] hbijm_old[k];
			delete[] hbijm_old;
		}
	}
	return;
}

void VMMC_CPUBackend::init() {
	MC_CPUBackend::init();

	// fix maxclust if evidently wrong
	if(_maxclust < 1) {
		OX_LOG(Logger::LOG_WARNING, "maxclust < 0, setting it to N = %i", N());
		_maxclust = N();
	}
	if(_maxclust > N()) {
		OX_LOG(Logger::LOG_WARNING, "maxclust > N does not make sense, setting it to N = %i", N());
		_maxclust = N();
	}

	if(_have_us) {
		_op.init_from_file(_op_file, _particles, N());
		_w.init((const char *) _weights_file, &_op, _safe_weights, _default_weight);
		if(_reload_hist) {
			_h.init(_init_hist_file, &_op, _etemps, _netemps);
		}
		else {
			_h.init(&_op, _etemps, _netemps);
		}
		_h.set_simtemp(_T);
	}

	if(_delta[MC_MOVE_TRANSLATION] * sqrt(3) > _verlet_skin) {
		throw oxDNAException("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");
	}

	_vmmc_box_side = _box->box_sides()[0];

	// setting the maximum displacement
	if(_preserve_topology) {
		_max_move_size = (number) 0.5;
		_max_move_size_sqr = _max_move_size * _max_move_size;
		OX_LOG(Logger::LOG_INFO, "Preserving topology; max_move_size = %lf...", _max_move_size);
	}
	else {
		_max_move_size = _vmmc_box_side / 2. - 2. * _rcut - 0.2;
		_max_move_size_sqr = _max_move_size * _max_move_size;
		OX_LOG(Logger::LOG_INFO, "Not attempting to preserve topology; max_move_size = %g", _max_move_size);
	}

	new_en3s.resize(N(), 0.);
	new_en5s.resize(N(), 0.);
	new_stn3s.resize(N(), 0.);
	new_stn5s.resize(N(), 0.);
	number tmpf, epq;
	BaseParticle * p, *q;
	for(int k = 0; k < N(); k++) {
		p = _particles[k];
		if(p->n3 != P_VIRTUAL) {
			q = p->n3;
			epq = _particle_particle_bonded_interaction_n3_VMMC(p, q, &tmpf);
			p->en3 = epq;
			q->en5 = epq;
			p->esn3 = tmpf;
			q->esn5 = tmpf;
		}
	}

	if(_small_system) {
		eijm = new number*[N()];
		eijm_old = new number*[N()];
		hbijm = new bool*[N()];
		hbijm_old = new bool*[N()];
		for(int k = 0; k < N(); k++) {
			eijm[k] = new number[N()];
			eijm_old[k] = new number[N()];
			hbijm[k] = new bool[N()];
			hbijm_old[k] = new bool[N()];
		}

		for(int k = 0; k < N(); k++) {
			for(int l = 0; l < k; l++) {
				p = _particles[k];
				q = _particles[l];
				if(p->n3 != q && p->n5 != q) {
					eijm[k][l] = eijm[l][k] = eijm_old[k][l] = eijm_old[l][k] = _particle_particle_nonbonded_interaction_VMMC(p, q, &tmpf);
					hbijm[k][l] = hbijm[l][k] = hbijm_old[k][l] = hbijm_old[l][k] = (tmpf < HB_CUTOFF);
				}
			}
		}
	}

	_init_cells();

	_compute_energy();

	check_overlaps();

	if(_have_us) {
		_update_ops();
		check_ops();
	}
}

void VMMC_CPUBackend::get_settings(input_file & inp) {
	MC_CPUBackend::get_settings(inp);
	int is_us, tmpi;

	CHECK_BOX("VMMC_CPUBackend", inp);

	std::string inter("");
	std::vector<std::string> ok_interactions; //= {"DNA", "DNA_nomesh", "DNA2", "DNA2_nomesh","RNA","RNA2"};//would work in C++11 and we woudln't to push back the elements
	ok_interactions.push_back("DNA");
	ok_interactions.push_back("DNA2");
	ok_interactions.push_back("DNA2ModInteraction");
	ok_interactions.push_back("DNA_nomesh");
	ok_interactions.push_back("DNA2_nomesh");
	ok_interactions.push_back("RNA");
	ok_interactions.push_back("RNA2");
	ok_interactions.push_back("NA");
	if(getInputString(&inp, "interaction_type", inter, 0) == KEY_FOUND) {
		// std::find points is equal to ok_interactions.end() if it can't find inter in ok_interactions.
		if(std::find(ok_interactions.begin(), ok_interactions.end(), inter) == ok_interactions.end()) {
			std::stringstream ok_as_string("");
			for(size_t i = 0; i < ok_interactions.size(); ++i)
				ok_as_string << ok_interactions[i] << ' ';
			throw oxDNAException("VMMC can be used only with the following interactions: %s", ok_as_string.str().c_str());
		}
	}

	if(getInputInt(&inp, "maxclust", &tmpi, 0) == KEY_FOUND) {
		_maxclust = tmpi;
		OX_LOG(Logger::LOG_INFO, "Using maxclust = %i", _maxclust);
	}

	if(getInputBoolAsInt(&inp, "small_system", &tmpi, 0) == KEY_FOUND) {
		if(tmpi > 0) {
			_small_system = true;
			OX_LOG(Logger::LOG_INFO, "Using algorithm N^2 for small system");
		}
		else {
			OX_LOG(Logger::LOG_INFO, "Using standard algorithm");
		}
	}

	getInputBool(&inp, "preserve_topology", &_preserve_topology, 0);

	if(getInputBoolAsInt(&inp, "umbrella_sampling", &is_us, 0) != KEY_NOT_FOUND) {
		if(is_us > 0) {
			_have_us = true;
			getInputString(&inp, "op_file", _op_file, 1);
			getInputString(&inp, "weights_file", _weights_file, 1);
			if(getInputString(&inp, "last_hist_file", _last_hist_file, 0) == KEY_NOT_FOUND) {
				sprintf(_last_hist_file, "last_hist.dat");
				OX_LOG(Logger::LOG_INFO, "Using default hist file %s", _last_hist_file);
			}
			if(getInputString(&inp, "traj_hist_file", _traj_hist_file, 0) == KEY_NOT_FOUND) {
				sprintf(_traj_hist_file, "traj_hist.dat");
				OX_LOG(Logger::LOG_INFO, "Using default traj hist file %s", _traj_hist_file);
			}

			// should we reload histograms?
			if(getInputString(&inp, "init_hist_file", _init_hist_file, 0) == KEY_FOUND) {
				OX_LOG(Logger::LOG_INFO, "Reloading histogram from %s", _init_hist_file);
				_reload_hist = true;
			}
			else {
				_reload_hist = false;
				FILE *temp_file = fopen(_traj_hist_file, "w");
				fclose(temp_file);
			}

			// whether to use unsafe weights
			if(getInputBoolAsInt(&inp, "safe_weights", &tmpi, 0) == KEY_FOUND) {
				_safe_weights = (tmpi > 0);
				if(!_safe_weights) {
					double tmpf;
					getInputDouble(&inp, "default_weight", &tmpf, 1);
					_default_weight = (number) tmpf;
					OX_LOG(Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) Using unsafe weights with default weight = %g", _default_weight);
				}
			}
			// whether to print out zero entries in traj_hist and last_hist
			if(getInputBoolAsInt(&inp, "skip_hist_zeros", &tmpi, 0) == KEY_FOUND) {
				_skip_hist_zeros = tmpi > 0;
				if(_skip_hist_zeros) {
					OX_LOG(Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) Skipping zero entries in traj_hist and last_hist files");
				}
			}

			// should we extrapolate the histogram at different
			// temperatures?
			char tstring[512];
			if(getInputString(&inp, "extrapolate_hist", tstring, 0) == KEY_FOUND) {
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
				}
				else {
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
		}
		else {
			_have_us = false;
		}
	}

	if(_have_us) {
		_h.read_interaction(inp);
	}

	if(getInputLLInt(&inp, "equilibration_steps", &_equilibration_steps, 0) == KEY_FOUND) {
		if(_equilibration_steps < 0) {
			OX_LOG(Logger::LOG_WARNING, "found equilibration_steps < 0; setting equilibration steps = 0");
			_equilibration_steps = 0;
		}
		OX_LOG(Logger::LOG_INFO, "(VMMC_CPUBackend.cpp) not recording histogram for steps < %lld", (long long unsigned)_equilibration_steps);
	}
}

// this function is just a wrapper that inverts p and q

inline number VMMC_CPUBackend::_particle_particle_bonded_interaction_n5_VMMC(BaseParticle *p, BaseParticle *q, number *stacking_en) {
	throw oxDNAException("ERROR: called a function that should not be called; file %s, line %d", __FILE__, __LINE__);
}

inline number VMMC_CPUBackend::_particle_particle_bonded_interaction_n3_VMMC(BaseParticle *p, BaseParticle *q, number *stacking_en) {
	BaseParticle * tmp1, *tmp2, *tmp3, *tmp4;
	tmp1 = p->n3;
	tmp2 = p->n5;
	tmp3 = q->n3;
	tmp4 = q->n5;
	p->n3 = q;
	p->n5 = P_VIRTUAL;
	q->n3 = P_VIRTUAL;
	q->n5 = p;

	if(stacking_en != 0)
	*stacking_en = 0;

	LR_vector r = q->pos - p->pos;

	// check for overlaps;
	LR_vector rback = r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
	number rbackr0;
	if(dynamic_cast<DNA2Interaction *>(_interaction.get()) != NULL) {
		rbackr0 = rback.module() - FENE_R0_OXDNA2;
	}
	else {
		rbackr0 = rback.module() - FENE_R0_OXDNA;
	}
	if(fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
		_overlap = true;
		p->n3 = tmp1;
		p->n5 = tmp2;
		q->n3 = tmp3;
		q->n5 = tmp4;
		return (number) (1.e6 * _T);
	}

	_interaction->set_computed_r(r);
	number energy = _interaction->pair_interaction_term(DNAInteraction::BACKBONE, p, q, false, false);
	energy += _interaction->pair_interaction_term(DNAInteraction::BONDED_EXCLUDED_VOLUME, p, q, false, false);
	number tmp_en = _interaction->pair_interaction_term(DNAInteraction::STACKING, p, q, false, false);
	energy += tmp_en;

	if(stacking_en != 0) {
		*stacking_en = tmp_en;
	}

	p->n3 = tmp1;
	p->n5 = tmp2;
	q->n3 = tmp3;
	q->n5 = tmp4;

	return energy;
}

inline number VMMC_CPUBackend::_particle_particle_nonbonded_interaction_VMMC(BaseParticle *p, BaseParticle *q, number *H_energy) {
	if(H_energy != 0)
	*H_energy = (number) 0;

	LR_vector r = _box->min_image(p->pos, q->pos);

	// early ejection
	if(r.norm() > _sqr_rcut) {
		return (number) 0.f;
	}

	_interaction->set_computed_r(r);
	number energy = _interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, p, q, false, false);

	if(H_energy != 0) {
		*H_energy = energy;
	}

	energy += _interaction->pair_interaction_term(DNAInteraction::NONBONDED_EXCLUDED_VOLUME, p, q, false, false);
	energy += _interaction->pair_interaction_term(DNAInteraction::CROSS_STACKING, p, q, false, false);

	// all interactions except DNA2Interaction use the DNAInteraction coaxial stacking*
	// *the hybrid interaction is a second exception
	if( (dynamic_cast<DNA2Interaction *>(_interaction.get()) == NULL) || (dynamic_cast<DRHInteraction *>(_interaction.get()) == NULL) ) {
		energy += _interaction->pair_interaction_term(DNAInteraction::COAXIAL_STACKING, p, q, false, false);
	}

	if(dynamic_cast<DRHInteraction *>(_interaction.get()) != NULL) {
		energy += _interaction->pair_interaction_term(DRHInteraction::COAXIAL_STACKING, p, q, false, false);
		energy += _interaction->pair_interaction_term(DRHInteraction::DEBYE_HUCKEL, p, q, false, false);
	}
	
	else if(dynamic_cast<DNA2Interaction *>(_interaction.get()) != NULL) {
		energy += _interaction->pair_interaction_term(DNA2Interaction::COAXIAL_STACKING, p, q, false, false);
		energy += _interaction->pair_interaction_term(DNA2Interaction::DEBYE_HUCKEL, p, q, false, false);
	}
	else if(dynamic_cast<RNA2Interaction *>(_interaction.get()) != NULL) {
		energy += _interaction->pair_interaction_term(RNA2Interaction::DEBYE_HUCKEL, p, q, false, false);
	}
	
	

	return energy;
}

inline bool find(int * clust, int size, int value) {
	int i;
	for(i = 0; i < size; i++) {
		if(clust[i] == value) {
			return true;
		}
	}
	return false;
}

inline void VMMC_CPUBackend::store_particle(BaseParticle * src) {
	BaseParticle *dst = _particles_old[src->index];

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->set_positions();
	dst->ext_potential = src->ext_potential;

	return;
}

inline void VMMC_CPUBackend::restore_particle(BaseParticle * dst) {
	BaseParticle *src = _particles_old[dst->index];

	dst->orientation = src->orientation;
	dst->orientationT = src->orientationT;
	dst->pos = src->pos;
	dst->set_positions();
	dst->ext_potential = src->ext_potential;

	return;
}

inline number VMMC_CPUBackend::build_cluster_small(movestr *moveptr, int maxsize, int *clust, int *size) {

	//printf ("before cycle...\n");
	//if (_have_us) check_ops();
	//printf ("passed\n");
	//printf ("\n\nfirst random: %g; %d @@@\n", drand48(), moveptr->seed);

	int nclust = 1;
	clust[0] = moveptr->seed;
	BaseParticle * pp, *qq;
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
	pp = _particles[clust[0]];
	pp->inclust = true;

	//ppold = &(_particles_old[pp->index]);
	store_particle(pp);
	_move_particle(moveptr, pp, pp);

	while(k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = _particles[clust[k]];

		//printf ("recruiting from %i... @@@\n", pp->index);

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if(_preserve_topology) {
			if(pp->pos.sqr_distance(_particles_old[pp->index]->pos) > _max_move_size_sqr) {
				_dU = 0.;
				_dU_stack = 0.;
				*size = nclust;
				_overlap = true;
				return 0.;
			}
		}

		// trying to recruit bonded neighbors of pp
		if(pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
			if(!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC(pp, qq, &stack_temp);
				test1 = VMMC_link(E_pp_moved, E_old);
				
				if(test1 > _next_rand()) {
					// prelink successful
					store_particle(qq);
					_move_particle(moveptr, qq, pp);

					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC(_particles_old[pp->index], qq);
					test2 = VMMC_link(E_qq_moved, E_old);
					if((test2 / test1) > _next_rand()) {
						//we did full_link, qq goes in the cluster
						//printf ("recruited %d ... @@@\n", qq->index);
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					}
					else {
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
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
		if(pp->n5 != P_VIRTUAL) {
			qq = pp->n5;
			if(!qq->inclust) {

				E_old = pp->en5;
				//E_pp_moved = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &stack_temp);
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC(qq, pp, &stack_temp);

				test1 = VMMC_link(E_pp_moved, E_old);
				if(test1 > _next_rand()) {
					// prelink successful
					store_particle(qq);
					_move_particle(moveptr, qq, pp);

					//E_qq_moved = _particle_particle_bonded_interaction_n5_VMMC (_particles_old[pp->index], qq);
					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC(qq, _particles_old[pp->index]);

					test2 = VMMC_link(E_qq_moved, E_old);
					if((test2 / test1) > _next_rand()) {
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

		number tmpf = (number) 0.;
		for(neigh = 0; neigh < N(); neigh++) {
			qq = _particles[neigh]; //qq is my neighbor
			if((!qq->inclust) && (pp->n3 != qq) && (pp->n5 != qq) && (qq != pp)) {

				E_old = eijm_old[pp->index][qq->index];

				if(E_old == (number) 0.)
				continue;

				E_pp_moved = _particle_particle_nonbonded_interaction_VMMC(pp, qq, &tmpf);
				test1 = VMMC_link(E_pp_moved, E_old);
				if(test1 > _next_rand()) {
					store_particle(qq);
					_move_particle(moveptr, qq, pp);

					E_qq_moved = _particle_particle_nonbonded_interaction_VMMC(_particles_old[pp->index], qq);

					test2 = VMMC_link(E_qq_moved, E_old);
					if((test2 / test1) > _next_rand()) {
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
						//printf ("recruited %d ... @@@\n", qq->index);
					}
					else {
						// prelinked;
						prelinked_particles.insert(qq->index);
						//_r_move_particle (moveptr, qq);
						restore_particle(qq);
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
		k++;
	}
	*size = nclust;

	if(nclust > maxsize) {
		pprime = 0.;
		_dU = 0.;
		_dU_stack = 0.;
		_overlap = true;
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
	for(int i = 0; i < nclust; i++) {
		prelinked_particles.erase(clust[i]);
	}
	if(prelinked_particles.size() > 0) {
		_reject_prelinks = true;
		_dU = 0;
		_dU_stack = 0;
		_overlap = true;
		//printf ("rejecting because of prelinks @@@");
		return (number) 0.;
	}

	// fix cells...
	for(int i = 0; i < nclust; i++) {
		int old_index, new_index;
		pp = _particles[clust[i]];
		old_index = _cells[pp->index];
		new_index = _get_cell_index(pp->pos);
		if(new_index != old_index) {
			_fix_list(pp->index, old_index, new_index);
		}
	}

	delta_E = delta_Est = E_anomaly = 0.;
	number tmpf_new, epq_new, epq_old;
	bool h_new, h_old;
	for(int i = 0; i < nclust; i++) {
		pp = _particles[clust[i]];
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

		for(neigh = 0; neigh < N(); neigh++) {
			qq = _particles[neigh]; //qq is my neighbor
			if((!qq->inclust) && (pp->n3 != qq) && (pp->n5 != qq) && (qq != pp)) {
				epq_old = eijm_old[pp->index][qq->index];

				if(epq_old == (number) 0.) {
					// in this case we have not considered it in the
					// recruiting stage...
					epq_new = _particle_particle_nonbonded_interaction_VMMC(pp, qq, &tmpf_new);
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
				if(epq_old == 0. && epq_new > 0.) {
					// we have just created an overlap where there
					// was no interaction
					E_anomaly -= epq_new;
				}

				// check for anomaly of first kind
				if(epq_old > 0. && epq_new == 0.) {
					// we have removed an overlap
					E_anomaly += epq_old;
				}

				// fix h_bonding...
				if(_have_us) {
					h_old = hbijm_old[pp->index][qq->index];
					h_new = hbijm[pp->index][qq->index];
					//if ((pp->index == 6 && qq->index == 9) || (pp->index == 9 && qq->index == 6)) printf ("ciao.. %i %i\n", h_old, h_new);
					if(h_old != h_new) {
						if(h_old == false) {
							_op.add_hb(pp->index, qq->index);
						}
						else {
							_op.remove_hb(pp->index, qq->index);
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

	pprime *= exp((1. / _T) * E_anomaly);

	//printf ("Returning %g @@@\n", pprime);

	_dU = delta_E;
	_dU_stack = delta_Est;

	return pprime;
}

// this function is the heart of the VMMC algorithm;this version uses
// cells to avoid the computation of O(N) interactions.

inline number VMMC_CPUBackend::build_cluster_cells(movestr *moveptr, int maxsize, int *clust, int *size) {
	int nclust = 1;
	clust[0] = moveptr->seed;
	BaseParticle * pp, *qq;
	number test1, test2;

	number pprime = 1.;

	number delta_E = 0;
	number delta_Est = 0;

	//check_ops();

	//_compute_energy();
	//_check_metainfo();
	//_check_old_metainfo();

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
	pp = _particles[clust[0]];
	pp->inclust = true;

	store_particle(pp);
	_move_particle(moveptr, pp, pp);

	assert(_overlap == false);

	// how far away am I recruiting?
	LR_vector how_far(0., 0., 0.);

	while(k < nclust && nclust <= maxsize) {
		//pp is the moved particle which is already in the cluster
		pp = _particles[clust[k]];

		// rejecting potential topology-breaking moves if appropiate key
		// is set in the input file
		if(_preserve_topology) {
			if(pp->pos.sqr_distance(_particles_old[pp->index]->pos) > _max_move_size_sqr) {
				_dU = 0.;
				_dU_stack = 0.;
				*size = nclust;
				_overlap = true;
				return 0.;
			}
		}

		// trying to recruit bonded neighbors of pp
		if(pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
			if(!qq->inclust) {
				//doing the prelinking
				E_old = pp->en3;
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC(pp, qq, &stack_temp);
				test1 = VMMC_link(E_pp_moved, E_old);
				if(_overlap || test1 > _next_rand()) {
					// prelink successful
					store_particle(qq);
					_move_particle(moveptr, qq, pp);

					// in case E_pp_moved created an overlap
					_overlap = false;

					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC(_particles_old[pp->index], qq);
					//_move_particle(moveptr, pp);

					test2 = VMMC_link(E_qq_moved, E_old);
					if(_overlap || (test2 / test1) > _next_rand()) {
						//we did full_link, qq goes in the cluster

						// in case E_qq_moved created an overlap
						_overlap = false;

						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					}
					else {
						assert(_overlap == false);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
						prelinked_particles.insert(qq->index);
					}
				}
				else {
					assert(_overlap == false);
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
		if(pp->n5 != P_VIRTUAL) {
			qq = pp->n5;
			if(!qq->inclust) {

				E_old = pp->en5;
				//E_pp_moved = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &stack_temp);
				E_pp_moved = _particle_particle_bonded_interaction_n3_VMMC(qq, pp, &stack_temp);

				test1 = VMMC_link(E_pp_moved, E_old);
				if(_overlap || test1 > _next_rand()) {
					// prelink successful
					store_particle(qq);
					_move_particle(moveptr, qq, pp);

					// in case we have recruited because of an overlap
					_overlap = false;

					//_r_move_particle(moveptr, pp);
					//E_qq_moved = _particle_particle_bonded_interaction_n5_VMMC (pp, qq);
					//_move_particle(moveptr, pp);
					//E_qq_moved = _particle_particle_bonded_interaction_n5_VMMC (_particles_old[pp->index], qq);
					E_qq_moved = _particle_particle_bonded_interaction_n3_VMMC(qq, _particles_old[pp->index]);

					test2 = VMMC_link(E_qq_moved, E_old);
					if(_overlap || (test2 / test1) > _next_rand()) {
						//we did full_link, qq goes to cluster
						_overlap = false;
						clust[nclust] = qq->index;
						qq->inclust = true;
						nclust++;
					}
					else {
						assert(_overlap == false);
						prelinked_particles.insert(qq->index);
						//_r_move_particle(moveptr, qq);
						restore_particle(qq);
					}
				}
				else {
					assert(_overlap == false);
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
		for(int c = 0; c < 27; c++) {
			icell = _neighcells[_cells[pp->index]][c];
			//icell = cell_neighbours (_cells[pp->index], c);
			//assert (icell == _neighcells[_cells[pp->index]][c]);
			neigh = _vmmc_heads[icell];
			while(neigh != P_INVALID) {
				qq = _particles[neigh]; //qq is my neighbor

				if(pp->n3 == qq || pp->n5 == qq) {
					neigh = qq->next_particle;
					continue;
				}

				if(qq->inclust == false) {
					E_old = _particle_particle_nonbonded_interaction_VMMC(_particles_old[pp->index], qq, &H_temp);

					if(E_old == (number) 0.) {
						neigh = qq->next_particle;
						continue;
					}

					E_pp_moved = _particle_particle_nonbonded_interaction_VMMC(pp, qq);

					test1 = VMMC_link(E_pp_moved, E_old);
					if(test1 > _next_rand()) {
						store_particle(qq);
						_move_particle(moveptr, qq, pp);

						//_r_move_particle (moveptr, pp);
						//E_qq_moved = _particle_particle_nonbonded_interaction_VMMC (pp, qq);
						//_move_particle (moveptr, pp);
						E_qq_moved = _particle_particle_nonbonded_interaction_VMMC(_particles_old[pp->index], qq);

						test2 = VMMC_link(E_qq_moved, E_old);
						if((test2 / test1) > _next_rand()) {
							clust[nclust] = qq->index;
							qq->inclust = true;
							nclust++;
						}
						else {
							// prelinked;
							prelinked_particles.insert(qq->index);
							//_r_move_particle (moveptr, qq);
							restore_particle(qq);
						}
					}
					else {
						if(fabs(E_old) > 0.) {
							// we store the possible interaction to account for later
							store_particle(qq);
							prev_inter.insert((pp->index > qq->index) ? (base_pair(qq->index, pp->index)) : (base_pair(pp->index, qq->index)));
						}
					}
				}
				neigh = qq->next_particle;
			}
		}
		k++;
	}
	*size = nclust;

	if(nclust > maxsize) {
		pprime = 0.;
		_dU = 0.;
		_dU_stack = 0.;
		_overlap = true;
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
	for(int i = 0; i < nclust; i++) {
		prelinked_particles.erase(clust[i]);
	}
	if(prelinked_particles.size() > 0) {
		//printf ("## setting pprime = 0. because of prelinked particles..\n");
		_reject_prelinks = true;
		_dU = 0;
		_dU_stack = 0;
		_overlap = true;
		return (number) 0.;
	}

	// fix cells...
	for(int i = 0; i < nclust; i++) {
		int old_index, new_index;
		pp = _particles[clust[i]];
		old_index = _cells[pp->index];
		new_index = _get_cell_index(pp->pos);
		if(new_index != old_index) {
			_fix_list(pp->index, old_index, new_index);
		}
	}

	delta_E = delta_Est = E_anomaly = 0.;
	number tmpf_old, tmpf_new, epq_new, epq_old;
	bool h_new, h_old;
	for(int i = 0; i < nclust; i++) {
		pp = _particles[clust[i]];
		if(pp->n3 != P_VIRTUAL) {
			qq = pp->n3;
			if(qq->inclust == false) {
				epq_new = _particle_particle_bonded_interaction_n3_VMMC(pp, qq, &tmpf_new);

				assert(_overlap == false);
				delta_E += epq_new - pp->en3;
				delta_Est += tmpf_new - pp->esn3;

				new_en3s[pp->index] = new_en5s[qq->index] = epq_new;
				new_stn3s[pp->index] = new_stn5s[qq->index] = tmpf_new;
			}
		}

		if(pp->n5 != P_VIRTUAL) {
			qq = pp->n5;
			if(qq->inclust == false) {
				if(_overlap)
				printf("cera prima\n");
				epq_new = _particle_particle_bonded_interaction_n3_VMMC(qq, pp, &tmpf_new);
				//epq_new = _particle_particle_bonded_interaction_n5_VMMC (pp, qq, &tmpf_new);

				assert(_overlap == false);

				delta_E += epq_new - pp->en5;
				delta_Est += tmpf_new - pp->esn5;

				new_en5s[pp->index] = new_en3s[qq->index] = epq_new;
				new_stn5s[pp->index] = new_stn3s[qq->index] = tmpf_new;
			}
		}

		for(int c = 0; c < 27; c++) {
			icell = _neighcells[_cells[pp->index]][c];
			//icell = cell_neighbours (_cells[pp->index], c);
			//assert (icell == _neighcells[_cells[pp->index]][c]);
			neigh = _vmmc_heads[icell];
			while(neigh != P_INVALID) {
				qq = _particles[neigh];

				if(pp->n3 == qq || pp->n5 == qq) {
					neigh = qq->next_particle;
					continue;
				}

				if(qq->inclust == false) {
					//_r_move_particle (moveptr, pp);
					//epq_old = _particle_particle_nonbonded_interaction_VMMC (pp, qq, &tmpf_old);
					//_move_particle (moveptr, pp);
					epq_old = _particle_particle_nonbonded_interaction_VMMC(_particles_old[pp->index], qq, &tmpf_old);
					epq_new = _particle_particle_nonbonded_interaction_VMMC(pp, qq, &tmpf_new);

					delta_E += epq_new - epq_old;

					// we have considered this interaction, so we remove it from the list
					if(fabs(epq_old) > 0.)
					prev_inter.erase((pp->index > qq->index) ? (base_pair(qq->index, pp->index)) : (base_pair(pp->index, qq->index)));

					// check for anomaly of second kind;
					if(epq_old == 0. && epq_new > 0.) {
						// we have just created an overlap where there
						// was no interaction
						E_anomaly -= epq_new;
					}

					// check for anomaly of first kind
					if(epq_old > 0.) {
						if(epq_new == 0.) {
							E_anomaly += epq_old;
						}
					}
					// fix h_bonding...
					if(_have_us) {
						h_new = tmpf_new < HB_CUTOFF;
						h_old = tmpf_old < HB_CUTOFF;
						//poss_breaks.erase ((pp->index>qq->index)?(base_pair(qq->index, pp->index)):(base_pair(pp->index, qq->index)));
						if(h_old != h_new) {
							if(h_old == false) {
								_op.add_hb(pp->index, qq->index);
							}
							else {
								_op.remove_hb(pp->index, qq->index);
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
	for(it = prev_inter.begin(); it != prev_inter.end(); it++) {
		number tmpf;
		pp = _particles[(*it).first];
		qq = _particles[(*it).second];
		if(!(pp->inclust && qq->inclust)) {
			epq_old = _particle_particle_nonbonded_interaction_VMMC(_particles_old[pp->index], _particles_old[qq->index], &tmpf);
			delta_E -= epq_old;
			if(epq_old > 0.)
			E_anomaly += epq_old;
			// if we get to this stage, epq_new has to be 0, so the hydrogen bond
			// is for sure not present anymore
			if(_have_us) {
				if(tmpf < HB_CUTOFF) {
					_op.remove_hb(pp->index, qq->index);
				}
			}
		}
	}

	pprime *= exp((1. / _T) * E_anomaly);

	_dU = delta_E;
	_dU_stack = delta_Est;

	return pprime;
}

inline void VMMC_CPUBackend::_move_particle(movestr *moveptr, BaseParticle *q, BaseParticle *p) {
	if(moveptr->type == MC_MOVE_TRANSLATION) {
		q->pos += moveptr->t;
	}
	else if(moveptr->type == MC_MOVE_ROTATION) {
		// in this case, the translation vector is the point around which we rotate
		BaseParticle *p_old = _particles_old[p->index];
		LR_vector dr, drp;
		if(p->strand_id == q->strand_id) {
			dr = q->pos - p_old->pos;
		}
		else {
			dr = _box->min_image(p_old->pos, q->pos);
		}

		drp = moveptr->R * dr;
		q->pos = p->pos + drp;
		q->orientation = moveptr->R * q->orientation;
		q->orientationT = q->orientation.get_transpose();
		q->set_positions();
	}
}

inline void VMMC_CPUBackend::_fix_list(int p_index, int oldcell, int newcell) {
	int j, jold;

	// remove p_index from its old cell
	//printf ("## %i %i %i\n", oldcell, _vmmc_N_cells, _cells[p_index]);
	j = _vmmc_heads[oldcell];
	jold = P_INVALID;
	assert(j != P_INVALID);
	while(j != p_index) {
		jold = j;
		j = _particles[j]->next_particle;
	}

	assert(j != jold);
	assert(j != P_INVALID);

	if(jold != P_INVALID) {
		_particles[jold]->next_particle = _particles[p_index]->next_particle;
	}
	else {
		_vmmc_heads[oldcell] = _particles[p_index]->next_particle;
	}

	// the old was full, since it contained the particle; is it now empty?
	if(_vmmc_heads[oldcell] == P_INVALID) {
		//printf ("now empty %i\n", oldcell);
		delete[] _neighcells[oldcell];
	}
	if(_vmmc_heads[newcell] == P_INVALID) {
		//printf ("now filled %i\n", newcell);
		_neighcells[newcell] = new int[27];
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
					nneigh++;
				}
			}
		}
	}

	// add it to the new cell
	_particles[p_index]->next_particle = _vmmc_heads[newcell];
	_vmmc_heads[newcell] = p_index;

	_cells[p_index] = newcell;

	return;
}

void VMMC_CPUBackend::sim_step() {

	_mytimer->resume();
	_timer_move->resume();

	LR_vector tmp;

	int *clust, nclust;
	clust = new int[N()];

	double oldweight, weight;
	int windex, oldwindex;
	oldweight = weight = 1.;
	if(_have_us) {
		oldweight = _w.get_weight(_op.get_all_states(), &oldwindex);
	}

	// set the potential due to external forces
	_U_ext = (number) 0.f;
	for(int k = 0; k < N(); k++) {
		BaseParticle *p = _particles[k];
		p->set_ext_potential(current_step(), _box.get());
		_U_ext += p->ext_potential;
	}

	for(int i = 0; i < N(); i++) {
		if(_have_us) {
			_op.store();
		}
		_dU_stack = 0.;

		// seed particle;
		int pi = (int) (drand48() * N());
		BaseParticle *p = _particles[pi];

		// select the move
		movestr move;
		move.seed = pi;
		move.type = (drand48() < 0.5) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;

		// generate translation / rotation
		if(move.type == MC_MOVE_TRANSLATION) {
			move.t = LR_vector(Utils::gaussian(), Utils::gaussian(), Utils::gaussian()) * _delta[MC_MOVE_TRANSLATION];
			move.R = LR_matrix((number) 1., (number) 0., (number) 0., (number) 0., (number) 1., (number) 0., (number) 0., (number) 0., (number) 1.);
		}
		else {
			//translation vector is then interpreted as the axis around
			//which we rotate by move_particle() below
			//pp = &(_particles[clust[0]]);
			move.R = Utils::get_random_rotation_matrix_from_angle(_delta[MC_MOVE_ROTATION] * Utils::gaussian());
			move.Rt = (move.R).get_transpose();
			//move.t = _particles[move.seed]->int_centers[DNANucleotide::BACK] + _particles[move.seed]->pos;
			move.t = _particles[move.seed]->int_centers[DNANucleotide::BACK];
		}
		_last_move = move.type;

		// build the cluster;
		//number pprime = build_cluster(pi, _maxclust, clust, &nclust, tainted, &ntainted);
		//printf("building cluster starting from %i...\n", move.seed);
		number pprime;
		if(_small_system) {
			pprime = build_cluster_small(&move, _maxclust, clust, &nclust);
		}
		else {
			pprime = build_cluster_cells(&move, _maxclust, clust, &nclust);
		}

		assert(nclust >= 1);

		number delta_E_ext = 0.;

		// if we are not SURE to reject the move, we check the external
		// forces. Otherwise, there is no point.
		if(_overlap == false && pprime > 0.) {
			for(int l = 0; l < nclust; l++) {
				p = _particles[clust[l]];
				delta_E_ext += -p->ext_potential;
				p->set_ext_potential(current_step(), _box.get());
				delta_E_ext += +p->ext_potential;
			}
			pprime *= exp(-(1. / _T) * delta_E_ext);
		}

		_op.fill_distance_parameters(_particles, _box.get());

		windex = oldwindex;
		weight = oldweight;
		if(_have_us) {
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
		//printf ("### %lf\n", _dU_stack);
		_tries[_last_move]++;

		// printf("## U: %lf dU: %lf, p': %lf, nclust: %d \n", _U, _dU, pprime, nclust);
		if(_overlap == false && pprime > drand48()) {
			if(nclust <= _maxclust)
			_accepted[_last_move]++;
			//if (!_reject_prelinks) _accepted[0]++;
			_U += _dU;
			_U_stack += _dU_stack;
			_U_ext += delta_E_ext;

			oldweight = weight; // if (!_have_us) oldweight = weight = 1.;
			oldwindex = windex; // if (!_have_us) oldweight = weight = 1.;

			for(int l = 0; l < nclust; l++) {
				BaseParticle * pp, *qq;
				pp = _particles[clust[l]];
				if(pp->n3 != P_VIRTUAL) {
					qq = pp->n3;
					if(qq->inclust == false) {
						pp->en3 = qq->en5 = new_en3s[clust[l]];
						pp->esn3 = qq->esn5 = new_stn3s[clust[l]];
					}
				}
				if(pp->n5 != P_VIRTUAL) {
					qq = pp->n5;
					if(qq->inclust == false) {
						pp->en5 = qq->en3 = new_en5s[clust[l]];
						pp->esn5 = qq->esn3 = new_stn5s[clust[l]];
					}
				}

				if(_small_system) {
					for(int c = 0; c < N(); c++) {
						qq = _particles[c];
						if(pp != qq && pp->n3 != qq && pp->n5 != qq && qq->inclust == false) {
							eijm_old[pp->index][qq->index] = eijm_old[qq->index][pp->index] = eijm[qq->index][pp->index];
							hbijm_old[pp->index][qq->index] = hbijm_old[qq->index][pp->index] = hbijm[pp->index][qq->index];
						}
					}
				}
			}
			//printf("## accepting dU = %lf, pprime = %lf\n", _dU, pprime);
			//printf("## checking metainfo after accepting\n");
			//_check_metainfo();
			//if (_have_us) check_ops();
			//printf("## checking metainfo after accepting: PASSED\n");
			//_check_old_metainfo();
		}
		else {
			//move rejected
			//printf("## rejecting dU = %lf, pprime = %lf, if %i==%i just updated lists\n", _dU, pprime, _just_updated_lists, true);
			for(int l = 0; l < nclust; l++) {
				int old_index, new_index;
				BaseParticle * pp;
				pp = _particles[clust[l]];
				//_r_move_particle (&move, pp);
				restore_particle(pp);
				old_index = _cells[pp->index];
				new_index = _get_cell_index(pp->pos);
				if(new_index != old_index) {
					_fix_list(pp->index, old_index, new_index);
				}
				pp->set_ext_potential(current_step(), _box.get());
			}

			_overlap = false;

			if(_have_us) {
				_op.restore();
			}
		}

		/*
		 // check ext potential
		 number c_ext = 0.;
		 number c_ext_fs = 0.;
		 for (int l = 0; l < N(); l ++) {
		 BaseParticle * pp;
		 pp = &(_particles[l]);
		 c_ext += pp->ext_potential;
		 pp->set_ext_potential(curr_step);
		 c_ext_fs += pp->ext_potential;
		 }
		 if (fabs (c_ext - c_ext_fs) > 1.e-6) {
		 fprintf (stderr, "%g %g -- beh\n", c_ext, c_ext_fs);
		 }
		 // check ext potential done*/

		// add to the histogram
		if(_have_us && current_step() > _equilibration_steps) {
			_h.add(oldwindex, oldweight, _U, _U_stack, _U_ext);
		}

		// reset the inclust property to the particles
		for(int k = 0; k < nclust; k++) {
			_particles[clust[k]]->inclust = false;
		}
	}

	//check_ops();

	delete[] clust;

	// check energy for percolation
	if(current_step() % (llint) _check_energy_every == 1) {
		//printf ("checking energy for percolation..\n");
		number U_from_tally = _U;
		_compute_energy();
		if((_U - U_from_tally) > 1.e-4) {
			throw oxDNAException("(VMMC_CPUBackend) Accumulated Energy (%g) and Energy computed from scratch (%g) don't match. Possibly percolating clusters. Your box is too small", U_from_tally, _U);
		}
		//printf ("all ok (%g %g)... \n", U_from_tally, _U);
	}

	_timer_move->pause();
	_mytimer->pause();
}

void VMMC_CPUBackend::check_ops() {
	if(!_have_us) {
		return;
	}

	int * state;

	//state = (int *) malloc(_op.get_hb_parameters_count() * sizeof(int));
	//memcpy(state, _op.get_hb_states(), _op.get_hb_parameters_count() * sizeof(int));
	state = (int *) malloc(_op.get_all_parameters_count() * sizeof(int));
	memcpy(state, _op.get_all_states(), _op.get_all_parameters_count() * sizeof(int));
	_op.reset();

	int i, j;
	BaseParticle *p, *q;
	number hpq;
	for(i = 0; i < N(); i++) {
		p = _particles[i];
		for(j = 0; j < i; j++) {
			q = _particles[j];
			if(p->n3 != q && p->n5 != q) {
				_particle_particle_nonbonded_interaction_VMMC(p, q, &hpq);
				if(hpq < HB_CUTOFF) {
					_op.add_hb(p->index, q->index);
				}
			}
		}
	}
	_op.fill_distance_parameters(_particles, _box.get());

	int * new_state = _op.get_all_states();
	int check = 0;

	for(i = 0; i < _op.get_all_parameters_count(); i++) {
		if(state[i] != new_state[i])
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

	if(check != 0) {
		printf("CASINO\n");
		abort();
	}

	free(state);
}

void VMMC_CPUBackend::_update_ops() {

	assert(_have_us);

	_op.reset();

	// hydrogen bonding
	int i, j, c;
	BaseParticle *p, *q;
	number hpq;
	for(i = 0; i < N(); i++) {
		p = _particles[i];
		for(c = 0; c < 27; c++) {
			j = _vmmc_heads[_neighcells[_cells[p->index]][c]];
			while(j != P_INVALID) {
				q = _particles[j];
				if(j < i && p->n3 != q && p->n5 != q) {
					_particle_particle_nonbonded_interaction_VMMC(p, q, &hpq);
					if(hpq < HB_CUTOFF) {
						_op.add_hb(i, j);
					}
				}
				j = q->next_particle;
			}
		}
	}

	// distances
	_op.fill_distance_parameters(_particles, _box.get());

	//exit(-1);
	return;
}

void VMMC_CPUBackend::_print_pos(int id) {
	BaseParticle * p;
	p = _particles[id];
	printf("%5i - np.array([% 10.6e % 10.6e % 10.6e])\n        np.array([% 10.6e % 10.6e % 10.6e])\n        np.array([% 10.6e % 10.6e % 10.6e])\n        np.array([% 10.6e % 10.6e % 10.6e])\n", id,
			p->pos.x, p->pos.y, p->pos.z, p->orientationT.v1.x, p->orientationT.v1.y, p->orientationT.v1.z, p->orientationT.v2.x, p->orientationT.v2.y, p->orientationT.v2.z, p->orientationT.v3.x,
			p->orientationT.v3.y, p->orientationT.v3.z);
	return;
}

inline void VMMC_CPUBackend::check_overlaps() {
	int i, N_overlaps = 0;
	//number epq;
	BaseParticle *p, *q;

	for(i = 0; i < N(); i++) {
		p = _particles[i];
		if(p->n3 != P_VIRTUAL) {
			_overlap = false;
			q = p->n3;
			_particle_particle_bonded_interaction_n3_VMMC(p, q);
			if(_overlap) {
				N_overlaps++;
				LR_vector rbb, r;
				r = p->pos - q->pos;
				rbb = r + p->int_centers[DNANucleotide::BACK] - q->int_centers[DNANucleotide::BACK];
				printf("### overlap %i and %i (%g, %g)\n", p->index, q->index, r.module(), rbb.module());
				_overlap = false;
			}
		}
	}
	if(N_overlaps > 0) {
		throw oxDNAException("VMMC: There is an overlap in the initial configuration.");
	}
}

char * VMMC_CPUBackend::get_op_state_str() {
	if(_have_us) {
		int * state = _op.get_all_states();
		char * aux;
		aux = (char *) _state_str;
		for(int i = 0; i < _op.get_all_parameters_count(); i++) {
			sprintf(aux, "%2d ", state[i]);
			aux = (char *) _state_str + strlen(_state_str);
		}
		sprintf(aux, "%lf", _w.get_weight(state));
		return _state_str;
	}
	else {
		sprintf(_state_str, " ");
		return _state_str;
	}
}

void VMMC_CPUBackend::print_conf(bool reduced, bool only_last) {
	SimBackend::print_conf(reduced, only_last);
	if(_have_us) {
		if(!only_last) {
			_h.print_to_file(_traj_hist_file, current_step(), false, _skip_hist_zeros);
		}
		_h.print_to_file(_last_hist_file, current_step(), true, _skip_hist_zeros);
	}
}

void VMMC_CPUBackend::print_conf(bool only_last) {
	SimBackend::print_conf(only_last);
	if(_have_us) {
		if(!only_last) {
			_h.print_to_file(_traj_hist_file, current_step(), false, _skip_hist_zeros);
		}
		_h.print_to_file(_last_hist_file, current_step(), true, _skip_hist_zeros);
	}
}

void VMMC_CPUBackend::_compute_energy() {
	_interaction->begin_energy_computation();
	// Since this function is called by MC_CPUBackend::init() but it uses cells initialized
	// by VMMC_CPUBackend::init(), we do nothing if it's called too early

	if(_vmmc_heads == NULL) {
		return;
	}
	BaseParticle * p, *q;
	_overlap = false;
	number res = (number) 0;
	number dres, tmpf;

	_U = (number) 0;

	for(int i = 0; i < N(); i++) {
		p = _particles[i];
		if(p->n3 != P_VIRTUAL) {
			q = p->n3;
			dres = _particle_particle_bonded_interaction_n3_VMMC(p, q);
			res += dres;
			_U += dres;
			if(_overlap) {
				printf("overlap found between particle %i and %i\n", p->index, q->index);
				_print_pos(2);
				_print_pos(1);
				abort();
			}
		}
		for(int c = 0; c < 27; c++) {
			int j = _vmmc_heads[_neighcells[_cells[p->index]][c]];
			while(j != P_INVALID) {
				q = _particles[j];
				if(p->n3 != q && p->n5 != q && p->index < q->index) {
					dres = _particle_particle_nonbonded_interaction_VMMC(p, q, &tmpf);
					_U += dres;
				}
				j = q->next_particle;
			}
		}
	}

	if(_overlap) {
		throw oxDNAException("overlap found. Aborting..\n");
	}
}

void VMMC_CPUBackend::_init_cells() {
	_vmmc_N_cells_side = (int) (floor(_vmmc_box_side / _rcut) + 0.1);
	while(_vmmc_N_cells_side > ceil(750.) && _vmmc_N_cells_side > 3) {
		_vmmc_N_cells_side--;
	}

	if(_vmmc_N_cells_side < 3)
	_vmmc_N_cells_side = 3;

	_vmmc_N_cells = _vmmc_N_cells_side * _vmmc_N_cells_side * _vmmc_N_cells_side;

	_vmmc_heads = new int[_vmmc_N_cells];
	_cells = new int[N()];
	_neighcells = new int *[_vmmc_N_cells];

	for(int i = 0; i < _vmmc_N_cells; i++)
		_vmmc_heads[i] = P_INVALID;

	for(int i = 0; i < N(); i++) {
		_particles[i]->next_particle = P_INVALID;
	}

	for(int i = 0; i < N(); i++) {
		BaseParticle *p = _particles[i];
		int cell_index = _get_cell_index(p->pos);
		int old_head = _vmmc_heads[cell_index];
		_vmmc_heads[cell_index] = i;
		_cells[i] = cell_index;
		assert(cell_index < _vmmc_N_cells);
		p->next_particle = old_head;
	}

	for(int i = 0; i < _vmmc_N_cells; i++) {
		if(_vmmc_heads[i] != P_INVALID) {
			_neighcells[i] = new int[27];
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
						nneigh++;
					}
				}
			}
		}
	}

	// check for the cells
	for(int i = 0; i < _vmmc_N_cells; i++) {
		//fprintf (stderr, "## i %i\n", i);
		int j = _vmmc_heads[i];
		//fprintf (stderr, "## j %i\n", j);
		if(j != P_INVALID) {
			// fprintf (stderr, "cell %5i: %i ", i, j);
			j = _particles[j]->next_particle;
			while(j != P_INVALID) {
				//fprintf (stderr, "%i ", j);
				j = _particles[j]->next_particle;
			}
			// fprintf (stderr, "\n");
		}
	}
	return;
}

void VMMC_CPUBackend::_delete_cells() {
	//for (int i = 0; i < _vmmc_N_cells; i ++) delete[] _neighcells[i];
	for(int i = 0; i < _vmmc_N_cells; i++) {
		if(_vmmc_heads[i] != P_INVALID) {
			delete[] _neighcells[i];
		}
	}
	delete[] _vmmc_heads;
	delete[] _cells;
	delete[] _neighcells;
	return;
}

void VMMC_CPUBackend::fix_diffusion() {
	// fix diffusion can sometimes change the value of the order paramer by changing the
	// orientations/coordinates particles that were barely above/below a cutoff.
	// To avoid this, we store the system before fix_diffusion and restore it 
	// if the order parameter value has changed.

	for(int i = 0; i < N(); i++)
		store_particle(_particles[i]);

	// check of order parameter
	std::vector<int> state(_op.get_all_parameters_count());
	state.insert(state.begin(), _op.get_all_states(), _op.get_all_states() + _op.get_all_parameters_count());
	_op.reset();

	SimBackend::fix_diffusion();

	_op.reset();
	for(int i = 0; i < N(); i++) {
		BaseParticle * p = _particles[i];
		for(int j = 0; j < i; j++) {
			BaseParticle * q = _particles[j];
			if(p->n3 != q && p->n5 != q) {
				number hpq;
				_particle_particle_nonbonded_interaction_VMMC(p, q, &hpq);
				if(hpq < HB_CUTOFF)
				_op.add_hb(p->index, q->index);
			}
		}
	}
	_op.fill_distance_parameters(_particles, _box.get());

	int *new_state = _op.get_all_states();

	int check = 0;
	for(int i = 0; i < _op.get_all_parameters_count(); i++) {
		check += abs(new_state[i] - state[i]);
	}

	if(check != 0) {
		OX_LOG(Logger::LOG_DEBUG, "(VMMC_CPUBackend) fix_diffusion() changed the value of the order parameter. Restoring simulation status before fix_diffusion()");
		for (int i = 0; i < N(); i ++) {
			restore_particle(_particles[i]);
		}
	}
}

void VMMC_CPUBackend::print_observables() {
	_lists->global_update(true);
	_backend_info += get_op_state_str();
	MCBackend::print_observables();
}

inline int VMMC_CPUBackend::_get_cell_index(const LR_vector &pos) {
	int res = (int) ((pos.x / _vmmc_box_side - floor(pos.x / _vmmc_box_side)) * (1. - std::numeric_limits<number>::epsilon()) * _vmmc_N_cells_side);
	res += _vmmc_N_cells_side * ((int) ((pos.y / _vmmc_box_side - floor(pos.y / _vmmc_box_side)) * (1. - std::numeric_limits<number>::epsilon()) * _vmmc_N_cells_side));
	res += _vmmc_N_cells_side * _vmmc_N_cells_side * ((int) ((pos.z / _vmmc_box_side - floor(pos.z / _vmmc_box_side)) * (1. - std::numeric_limits<number>::epsilon()) * _vmmc_N_cells_side));
	return res;
}

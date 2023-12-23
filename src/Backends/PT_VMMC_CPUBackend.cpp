/*
 * PT_VMMC_CPUBackend.cpp
 *
 *  Created on: 26/sep/2012
 *      Author: Flavio
 */

#include "PT_VMMC_CPUBackend.h"

#include <mpi.h>

PT_VMMC_CPUBackend::PT_VMMC_CPUBackend() :
				VMMC_CPUBackend() {
	// parallel tempering
	_npttemps = 4;
	_pttemps = NULL;
	_my_mpi_id = -1;
	_pt_move_every = 1000;
	_pt_exchange_tries = (llint) 1;
	_pt_exchange_accepted = (llint) 1;
	_pt_common_weights = false;
	_U_ext = (number) 0.;
}

PT_VMMC_CPUBackend::~PT_VMMC_CPUBackend() {
	delete[] _exchange_conf;
	delete[] _pttemps;
}

void PT_VMMC_CPUBackend::init() {
	if(_oxRNA_stacking) {
		RNAInteraction *it = dynamic_cast<RNAInteraction*>(_interaction.get());
		model = it->get_model();
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &(_my_mpi_id));
	MPI_Comm_size(MPI_COMM_WORLD, &(_mpi_nprocs));

	char my_conf_filename[1024];
	sprintf(my_conf_filename, "%s%d", _conf_filename.c_str(), _my_mpi_id);

	_conf_filename = string(my_conf_filename);

	// check that temperatures are in order...
	bool check2 = true;
	if(_my_mpi_id == 0) {
		for(int k = 0; k < _npttemps - 1; k++) {
			if(_pttemps[k] >= _pttemps[k + 1]) {
				check2 = false;
			}
		}
		if(!check2) {
			fprintf(stderr, "Attention: temperatures are not in increasing order. Hoping for the best\n");
			//OX_LOG(Logger::LOG_INFO, "Attention: temperatures are not in increasing order. Hoping for the best");
		}
	}

	// let's get our own temperature...
	if(_npttemps != _mpi_nprocs) {
		throw oxDNAException("Number of PT temperatures does not match number of processes (%d != %d)", _npttemps, _mpi_nprocs);
	}
	_T = _pttemps[_my_mpi_id];
	//OX_LOG(Logger::LOG_INFO, "Replica %d: Running at T=%g", _my_mpi_id, _T);
	fprintf(stderr, "Replica %d: Running at T=%g\n", _my_mpi_id, _T);

	CONFIG_INFO->update_temperature(_T);

	fprintf(stderr, "REPLICA %d: reading configuration from %s\n", _my_mpi_id, my_conf_filename);
	VMMC_CPUBackend::init();
	// N() returns no non-sense only after having called init()
	_exchange_conf = new PT_serialized_particle_info[N()];

	//fprintf (stderr, "REPLICA %d: Running at T=%g\n", _my_mpi_id, _T);

	// setting replica number
	_which_replica = _my_mpi_id;

	// changing filename
	char extra[16];
	sprintf(extra, "%d", _my_mpi_id);

	if(_reload_hist) {
		strcat(_init_hist_file, extra);
	}

	// common weights file? if so, we have a single file
	if(_have_us) {
		strcat(_last_hist_file, extra);
		strcat(_traj_hist_file, extra);

		sprintf(_irresp_weights_file, "%s", _weights_file);

		if(_pt_common_weights == false) {
			strcat(_weights_file, extra);

			if(_my_mpi_id != (_mpi_nprocs - 1)) {
				char extra2[16];
				sprintf(extra2, "%d", _my_mpi_id + 1);
				strcat(_irresp_weights_file, extra2);
			}
		}

		_irresp_w.init((const char*) _irresp_weights_file, &_op, _safe_weights, _default_weight);

		fprintf(stderr, "(from replica %d) common_weights = %d; weights file = %s\n", _my_mpi_id, _pt_common_weights, _weights_file);
	}
}

void PT_VMMC_CPUBackend::get_settings(input_file &inp) {
	VMMC_CPUBackend::get_settings(inp);
	// let's get the temperatures at which to run PT
	char tstring[512];

	getInputInt(&inp, "pt_every", &_pt_move_every, 0);

	if(getInputString(&inp, "pt_temp_list", tstring, 1) == KEY_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Running Parallel Tempering at temperatures .... %s\n", tstring);
		char *aux, deg;
		int c = 0, check;
		double *tmpt;
		tmpt = new double[100];
		aux = strtok(tstring, ",");
		while(aux != NULL) {
			//printf ("parsing %s\n", aux);
			check = sscanf(aux, "%lf %c", &(tmpt[c]), &deg);
			if(check < 1) {
				throw oxDNAException("Unrecognizable line in pt_temp_list");
			}
			if(check == 1) {
				; // do nothing
			}
			if(check == 2) {
				deg = tolower(deg);
				switch(deg) {
				case 'c':
					tmpt[c] = (tmpt[c] + 273.15) * 0.1 / 300.;
					break;
				case 'k':
					tmpt[c] = tmpt[c] * 0.1 / 300.;
					break;
				default:
					throw oxDNAException("Unrecognizable temperature '%s' in pt_temp_list", tmpt[c]);
					break;
				}
			}
			c++;
			aux = strtok(NULL, ",");
		}
		if(c == 0) {
			throw oxDNAException("Nothing found in pt_temp_list");
		}
		else {
			fprintf(stderr, "pt_temp_list: ");
			for(int i = 0; i < c; i++)
				fprintf(stderr, "%lf ", tmpt[i]);
			fprintf(stderr, "\n");
		}

		_npttemps = c;
		if(_npttemps > 0) {
			_pttemps = new double[_npttemps];
			memcpy(_pttemps, tmpt, _npttemps * sizeof(double));
		}
		delete[] tmpt;
	}

	int tmpi = -1;
	if(getInputInt(&inp, "pt_common_weights", &tmpi, 0) == KEY_FOUND) {
		if(tmpi <= 0)
			_pt_common_weights = false;
		else
			_pt_common_weights = true;
	}

	// For the stacking energy calculations
	_oxDNA2_stacking = false;
	_oxRNA_stacking = false;

	std::string inter_type("");
	getInputString(&inp, "interaction_type", inter_type, 0);
	if(inter_type.compare("DNA2") == 0) {
		_oxDNA2_stacking = true;
	}
	if(inter_type.compare("RNA2") == 0 || inter_type.compare("RNA2") == 0) {
		_oxRNA_stacking = true;
	}
}

void PT_serialized_particle_info::read_from(BaseParticle *par) {
	pos = par->pos;
	orientation = par->orientation;
}

void PT_serialized_particle_info::write_to(BaseParticle *par) {
	par->pos = pos;
	par->orientation = orientation;
}

void PT_VMMC_CPUBackend::sim_step() {
	//printf("This is a step %ld in process %d, with ene %f \n",curr_step,_my_mpi_id,_U);
	VMMC_CPUBackend::sim_step();

	if(current_step() % _pt_move_every == 0 && current_step() > 2) {

		// if we have forces, we should compute the external potential
		_U_ext = (number) 0.;
		if(_external_forces) {
			BaseParticle *p;
			for(int i = 0; i < N(); i++) {
				p = _particles[i];
				p->set_ext_potential(current_step(), _box.get());
				_U_ext += p->ext_potential;
			}
		}

		// find out if we try odd or even pairs
		//printf ("(from %d) attempting move exchange...\n", _my_mpi_id);
		bool odd_pairs = (((current_step() / _pt_move_every) % 2) == 0);
		bool im_responsible = ((_my_mpi_id % 2) == odd_pairs);
		int resp_id, irresp_id;
		if(im_responsible) {
			resp_id = _my_mpi_id;
			irresp_id = _my_mpi_id + 1;
			//printf ("(from %d); this turn (curr_step = %u) Im responsible of myself and %i\n", _my_mpi_id, (unsigned) curr_step, irresp_id);
		}
		else {
			resp_id = _my_mpi_id - 1;
			irresp_id = _my_mpi_id;
			//printf ("(from %d); this turn (curr_step = %u) Im NOT responsible; I rely on %i\n", _my_mpi_id, (unsigned) curr_step, resp_id);
		}

		// if I am responsible, I wait for the energy and conf of the next guy,
		// read it, decide whether to accept the move, send the conf back;
		if(im_responsible) {
			// check if the next guy exists:
			if(irresp_id >= 0 && irresp_id < _mpi_nprocs) {
				_pt_exchange_tries++;
				// Marvin Waitforit Eriksen
				_get_exchange_energy(irresp_id);
				//printf ("(from %d) got %f as energy...\n", _my_mpi_id, _exchange_energy.U);
				_get_exchange_conf(irresp_id);
				// make up my mind
				number fact, b1u1, b2u2, b2u1, b1u2, b1, b2, et, e0;

				// the fact that the potential is temperature dependent introduces this horrible
				// extrapolating procedure
				b1 = 1. / _T;
				b2 = 1. / _exchange_energy.T;

				if(_oxDNA2_stacking)
					et = _U_stack * STCK_FACT_EPS_OXDNA2 / (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 / b1);
				else if(_oxRNA_stacking)
					et = _U_stack * model->RNA_STCK_FACT_EPS / (model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS / b1);
				else
					et = _U_stack * STCK_FACT_EPS_OXDNA / (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA / b1);

				e0 = _U - et / b1;

				b2u1 = b2 * (e0 + et / b2);

				if(_oxDNA2_stacking)
					et = _exchange_energy.U_stack * STCK_FACT_EPS_OXDNA2 / (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 / b2);
				else if(_oxRNA_stacking)
					et = _exchange_energy.U_stack * model->RNA_STCK_FACT_EPS / (model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS / b2);
				else
					et = _exchange_energy.U_stack * STCK_FACT_EPS_OXDNA / (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA / b2);
				e0 = _exchange_energy.U - et / b2;

				b1u2 = b1 * (e0 + et / b1);

				b1u1 = b1 * _U;
				b2u2 = b2 * _exchange_energy.U;

				fact = exp(b1u1 + b2u2 - b2u1 - b1u2);

				// let's take the weights into account
				if(_have_us) {
					number w21, w12, w22, w11;

					w11 = _w.get_weight(_op.get_all_states());
					w22 = _exchange_energy.w;

					// weight of my conf in the irresponsible simulation
					w21 = _irresp_w.get_weight(_op.get_all_states());
					// weight of the other conf in this simulation:
					// we do something slightly obscure: we exchange the weight index
					// for convenience. Go in Weights.cpp if you want to find out/remember
					// what it does
					w12 = _w.get_weight_by_index(_exchange_energy.weight_index);

					fact *= ((w21 * w12) / (w22 * w11));
				}

				// we should take care of the external potential here...
				if(_external_forces) {
					fact *= exp((b1 - b2) * (_U_ext - _exchange_energy.U_ext));
				}

				//printf ("(from %d) fact: %lf\n", _my_mpi_id, fact);
				if(drand48() < fact) {
					_pt_exchange_accepted++;
					// send my conf over there
					// store the other guy's
					//printf ("(from %d) accepting exchange...\n", _my_mpi_id);
					PT_serialized_particle_info *buffer_conf;
					PT_energy_info buffer_energy;
					buffer_conf = new PT_serialized_particle_info[N()];
					for(int i = 0; i < N(); i++) {
						buffer_conf[i].pos = _exchange_conf[i].pos;
						buffer_conf[i].orientation = _exchange_conf[i].orientation;
					}
					buffer_energy.U = _exchange_energy.U;
					buffer_energy.U_stack = _exchange_energy.U_stack;
					buffer_energy.T = _exchange_energy.T;
					buffer_energy.U_ext = _exchange_energy.U_ext;
					buffer_energy.replica_id = _exchange_energy.replica_id;

					_build_exchange_energy();
					_send_exchange_energy(irresp_id);

					// put my current one in exchange_conf
					_build_exchange_conf();
					_send_exchange_conf(irresp_id);

					// put the other guy's old conf in my echange vector
					for(int i = 0; i < N(); i++) {
						_exchange_conf[i].pos = buffer_conf[i].pos;
						_exchange_conf[i].orientation = buffer_conf[i].orientation;
					}
					_exchange_energy.U = buffer_energy.U;
					_exchange_energy.U_stack = buffer_energy.U_stack;
					_exchange_energy.T = buffer_energy.T;
					_exchange_energy.U_ext = buffer_energy.U_ext;
					_exchange_energy.replica_id = buffer_energy.replica_id;

					// rebuild my conf (which will now become the other guy's old one
					_rebuild_exchange_conf();
					_rebuild_exchange_energy();

					delete[] buffer_conf;

				}
				else {
					//printf ("(from %d) rejecting exchange...\n", _my_mpi_id);
					// send back the guy's conf as it was
					_send_exchange_energy(irresp_id);
					_send_exchange_conf(irresp_id);
					// send back either its old one or my old one;
				}
				VMMC_CPUBackend::_compute_energy();
			}
		}
		else {
			// not responsible. I basically don't know anything, but perhaps
			// more importantly I have no interest in doing so.
			if(resp_id >= 0 && resp_id < _mpi_nprocs) {
				// send my energy and conf
				_build_exchange_energy();
				_send_exchange_energy(resp_id);

				_build_exchange_conf();
				_send_exchange_conf(resp_id);

				// wait for the conf to use next;
				_get_exchange_energy(resp_id);
				_get_exchange_conf(resp_id);

				_rebuild_exchange_energy();
				_rebuild_exchange_conf();
				VMMC_CPUBackend::_compute_energy();
				//printf("Received the following positions in process %d\n",_my_mpi_id);
				//for (int kk = 0; kk < N(); kk++)
				//	{
				//	_print_pos(kk);
				//}
				fflush(stdout);
			}
		}
		fflush(stdout);
		// debug
		if(fabs(_T - _pttemps[_my_mpi_id]) > 1e-7) {
			fprintf(stderr, "DISASTRO\n\n\n");
		}
		// we should set the forces again if we have swapped conf
		if(_external_forces) {
			BaseParticle *p;
			for(int i = 0; i < N(); i++) {
				p = _particles[i];
				p->set_ext_potential(current_step(), _box.get());
			}
		}

	}
}

void PT_VMMC_CPUBackend::_send_exchange_energy(int other_id) {
	//printf ("(from %d) sending energy info to %d\n", _my_mpi_id, other_id);
	_MPI_send_block_data((void*) (&_exchange_energy), sizeof(PT_energy_info), other_id);
}

void PT_VMMC_CPUBackend::_get_exchange_energy(int other_id) {
	//printf ("(from %d) waiting energy info from %d\n", _my_mpi_id, other_id);
	_MPI_receive_block_data((void*) (&_exchange_energy), sizeof(PT_energy_info), other_id);
}

void PT_VMMC_CPUBackend::_send_exchange_conf(int other_id) {
	//printf ("(from %d) sending conf info to %d\n", _my_mpi_id, other_id);
	_MPI_send_block_data((void*) (_exchange_conf), N() * sizeof(PT_serialized_particle_info), other_id);
}

void PT_VMMC_CPUBackend::_get_exchange_conf(int other_id) {
	//printf ("(from %d) waiting conf info from %d\n", _my_mpi_id, other_id);
	_MPI_receive_block_data((void*) (_exchange_conf), N() * sizeof(PT_serialized_particle_info), other_id);
}

void PT_VMMC_CPUBackend::_build_exchange_conf() {
	for(int i = 0; i < N(); i++) {
		_exchange_conf[i].read_from(_particles[i]);
	}
	return;
}

void PT_VMMC_CPUBackend::_build_exchange_energy() {
	_exchange_energy.U = _U;
	_exchange_energy.U_stack = _U_stack;
	_exchange_energy.T = _T;
	_exchange_energy.U_ext = _U_ext;
	_exchange_energy.replica_id = _which_replica;
	if(_have_us)
		_exchange_energy.w = _w.get_weight(_op.get_all_states(), &(_exchange_energy.weight_index));
}

void PT_VMMC_CPUBackend::_rebuild_exchange_energy() {
	if(_oxDNA2_stacking)
		_U_stack = _exchange_energy.U_stack * (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * _T) / (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * _exchange_energy.T);
	else if(_oxRNA_stacking)
		_U_stack = _exchange_energy.U_stack * (model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS * _T) / (model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS * _exchange_energy.T);
	else
		_U_stack = _exchange_energy.U_stack * (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * _T) / (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * _exchange_energy.T);

	_U = _exchange_energy.U + _U_stack - _exchange_energy.U_stack;
	_which_replica = _exchange_energy.replica_id;
	// N.B. we ignore the temperature, as we need it to decide wether to accept
	// or not, but we don't want to overwrite our current one
}

void PT_VMMC_CPUBackend::_rebuild_exchange_conf() {
	BaseParticle *p;
	for(int i = 0; i < N(); i++) {
		p = _particles[i];
		_exchange_conf[i].write_to(p);
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();
	}
	VMMC_CPUBackend::_delete_cells();
	VMMC_CPUBackend::_init_cells();

	number tmpf, epq;
	BaseParticle *q;
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

	// here we reset order parameters
	_op.reset();
	int i, j;
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

	//VMMC_CPUBackend::_update_metainfo();
	return;
}

int PT_VMMC_CPUBackend::_MPI_send_block_data(void *data, size_t size, int node_to, int TAG) {
	int ret = MPI_Send((void*) data, size, MPI_CHAR, node_to, TAG, MPI_COMM_WORLD);
	if(ret != MPI_SUCCESS) {
		throw oxDNAException("Error while sending MPI message");
		return -1;
	}
	else {
		return 1;
	}
}

int PT_VMMC_CPUBackend::_MPI_receive_block_data(void *data, size_t size, int node_from, int TAG) {
	MPI_Status stat;
	int ret = MPI_Recv((void*) data, size, MPI_CHAR, node_from, TAG, MPI_COMM_WORLD, &stat);

	if(ret != MPI_SUCCESS) {
		throw oxDNAException("Error while receving MPI message");
		return -1;
	}
	else {
		return 1;
	}
}

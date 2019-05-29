/*
 * PT_VMMC_CPUBackend.cpp
 *
 *  Created on: 26/sep/2012
 *      Author: Flavio
 */

#include "PT_VMMC_CPUBackend.h"

template<typename number> PT_VMMC_CPUBackend<number>::PT_VMMC_CPUBackend() : VMMC_CPUBackend<number> () {
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

template<typename number>
PT_VMMC_CPUBackend<number>::~PT_VMMC_CPUBackend() {
	delete [] _exchange_conf;
	delete [] _pttemps;
}

template<typename number>
void PT_VMMC_CPUBackend<number>::init () {
	//VMMC_CPUBackend<number>::init(conf_input);

	MPI_Comm_rank (MPI_COMM_WORLD, &(_my_mpi_id));
	MPI_Comm_size (MPI_COMM_WORLD, &(_mpi_nprocs));

	char my_conf_filename[256];
	sprintf (my_conf_filename, "%s%d", conf_filename, _my_mpi_id);

	fprintf (stderr, "REPLICA %d: reading configuration from %s\n", _my_mpi_id, my_conf_filename);
	VMMC_CPUBackend<number>::init (my_conf_filename);

	_exchange_conf = new PT_serialized_particle_info<number>[this->_N];

	// check that temperatures are in order...
	bool check2 = true;
	if (this->_my_mpi_id == 0) {
		for (int k = 0; k < _npttemps - 1; k++) {
			if (_pttemps[k] >= _pttemps[k + 1]) {
				check2 = false;
			}
		}
		if (!check2) {
			fprintf (stderr, "Attention: temperatures are not in increasing order. Hoping for the best\n");
			//OX_LOG(Logger::LOG_INFO, "Attention: temperatures are not in increasing order. Hoping for the best");
		}
	}

	// let's get our own temperature...
	if (_npttemps != _mpi_nprocs) throw oxDNAException("Number of PT temperatures does not match number of processes (%d != %d)", _npttemps, _mpi_nprocs);
	this->_T = _pttemps[_my_mpi_id];
	//OX_LOG(Logger::LOG_INFO, "Replica %d: Running at T=%g", _my_mpi_id, this->_T);

	OX_LOG(Logger::LOG_INFO, "Deleting previous interaction and creating one...");
	//this->_interaction.init(this->_T);
	throw oxDNAException ("File %s, line %d: not implemented", __FILE__, __LINE__);
	// here we should initialize the interaction to use the appropriate T
	
	OX_LOG(Logger::LOG_INFO, "Deleting previous lists and creating one...");
	throw oxDNAException ("File %s, line %d: not implemented", __FILE__, __LINE__);
	// here we should create our own lists...

	//VMMC_CPUBackend<number>::_compute_energy();
	MC_CPUBackend<number>::_compute_energy();

	fprintf (stderr, "REPLICA %d: Running at T=%g\n", _my_mpi_id, this->_T);

	// setting replica number
	_which_replica = _my_mpi_id;

	// changing filename
	char extra[16];
	sprintf (extra, "%d", _my_mpi_id);
	strcat (this->_last_hist_file, extra);
	if (this->_reload_hist) strcat (this->_init_hist_file, extra);

	// common weights file? if so, we have a single file
	if (this->_have_us) {

		sprintf (_irresp_weights_file, "%s", this->_weights_file);

		if (_pt_common_weights == false) {
			strcat (this->_weights_file, extra);

			if (_my_mpi_id != (_mpi_nprocs - 1)) {
				char extra2[16];
				sprintf (extra2, "%d", _my_mpi_id + 1);
				strcat (_irresp_weights_file, extra2);
			}
		}

		_irresp_w.init((const char *) _irresp_weights_file, &this->_op);

		fprintf (stderr, "(from replica %d) common_weights = %d; weights file = %s\n", _my_mpi_id, _pt_common_weights, this->_weights_file);
	}
}

template<typename number>
void PT_VMMC_CPUBackend<number>::get_settings(input_file & inp) {
	VMMC_CPUBackend<number>::get_settings(inp);
	// let's get the temperatures at which to run PT
	char tstring[512];

	getInputInt(&inp, "pt_every", &_pt_move_every, 0);

	if (getInputString(&inp, "pt_temp_list", tstring, 1) == KEY_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Running Parallel Tempering at temperatures .... %s\n", tstring);
		char * aux, deg;
		int c = 0, check;
		double * tmpt;
		tmpt = new double[100];
		aux = strtok(tstring, ",");
		while (aux != NULL) {
			//printf ("parsing %s\n", aux);
			check = sscanf(aux, "%lf %c", &(tmpt[c]), &deg);
			if (check < 1) {
				throw oxDNAException("Unrecognizable line in pt_temp_list");
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
					throw oxDNAException("Unrecognizable temperature '%s' in pt_temp_list", tmpt[c]);
					break;
				}
			}
			c++;
			aux = strtok(NULL, ",");
		}
		if (c == 0) {
			throw oxDNAException("Nothing found in pt_temp_list");
		} else {
			fprintf(stderr, "pt_temp_list: ");
			for (int i = 0; i < c; i++)
				fprintf(stderr, "%lf ", tmpt[i]);
			fprintf(stderr, "\n");
		}

		_npttemps = c;
		if (_npttemps > 0) {
			_pttemps = new double[_npttemps];
			memcpy(_pttemps, tmpt, _npttemps * sizeof(double));
		}
		delete[] tmpt;
		//abort ();
	}

	int tmpi = -1;
	if (getInputInt(&inp, "pt_common_weights", &tmpi, 0) == KEY_FOUND) {
		if (tmpi <= 0) _pt_common_weights = false;
		else _pt_common_weights = true;
	}

	// For the stacking energy calculations
	_oxDNA2_stacking = false;
	std::string inter_type ("");
	getInputString(&inp, "interaction_type", inter_type, 0);
	if(inter_type.compare("DNA2") == 0) _oxDNA2_stacking = true;
}

template <typename number>
void PT_serialized_particle_info<number>::read_from (BaseParticle<number> * par) {
	pos = par->pos;
	orientation= par->orientation;
}

template <typename number>
void PT_serialized_particle_info<number>::write_to (BaseParticle<number> * par) {
	par->pos = pos;
	par->orientation = orientation;
}

template<typename number>
void PT_VMMC_CPUBackend<number>::sim_step(llint curr_step) {
	VMMC_CPUBackend<number>::sim_step(curr_step);

	if (curr_step % _pt_move_every == 0 && curr_step > 2) {

		// if we have forces, we should compute the external potential
		_U_ext = (number) 0.;
		if (this->_external_forces) {
			BaseParticle<number> * p;
			for (int i = 0; i < this->_N; i ++) {
				p = this->_particles[i];
				p->set_ext_potential(curr_step, this->_box_side);
				_U_ext += p->ext_potential;
			}
		}

		// find out if we try odd or even pairs
		//printf ("(from %d) attempting move exchange...\n", _my_mpi_id);
		bool odd_pairs = (((curr_step / _pt_move_every) % 2) == 0);
		bool im_responsible = ((_my_mpi_id%2) == odd_pairs);
		int resp_id, irresp_id;
		if (im_responsible) {
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
		if (im_responsible ) {
			// check if the next guy exists:
			if (irresp_id >= 0 && irresp_id < _mpi_nprocs) {
				this->_pt_exchange_tries ++;
				// Marvin Waitforit Eriksen
				_get_exchange_energy (irresp_id);
				//printf ("(from %d) got %f as energy...\n", _my_mpi_id, _exchange_energy.U);
				_get_exchange_conf (irresp_id);
				// make up my mind
				number fact, b1u1, b2u2, b2u1, b1u2, b1, b2, et, e0;

				// the fact that the potential is temperature dependent introduces this horrible
				// extrapolating procedure
				b1 = 1./this->_T;
				b2 = 1./_exchange_energy.T;

				if (_oxDNA2_stacking) et = this->_U_stack * STCK_FACT_EPS_OXDNA2 / (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 / b1);
				else et = this->_U_stack * STCK_FACT_EPS_OXDNA / (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA / b1);

				e0 = this->_U - et / b1;

				b2u1 = b2 * (e0 + et / b2);

				if (_oxDNA2_stacking) et = _exchange_energy.U_stack * STCK_FACT_EPS_OXDNA2 / (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 / b2);
				else et = _exchange_energy.U_stack * STCK_FACT_EPS_OXDNA / (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA / b2);
				e0 = _exchange_energy.U - et / b2;

				b1u2 = b1 * (e0 + et / b1);

				b1u1 = b1 * this->_U;
				b2u2 = b2 * _exchange_energy.U;

				fact = exp (b1u1 + b2u2 - b2u1 - b1u2);

				// let's take the weights into account
				if (this->_have_us) {
					number w21, w12, w22, w11;

					w11 = this->_w.get_weight (this->_op.get_hb_states());
					w22 = this->_exchange_energy.w;

					// weight of my conf in the irresponsible simulation
					w21 = _irresp_w.get_weight (this->_op.get_hb_states());
					// weight of the other conf in this simulation:
					// we do something slightly obscure: we exchange the weight index
					// for convenience. Go in Weights.cpp if you want to find out/remember
					// what it does
					w12 = this->_w.get_weight_by_index (_exchange_energy.weight_index);

					fact *= ((w21 * w12) / (w22 * w11));
				}
				//fprintf (stderr, "(replica %d) had %lf %lf %lf AAA\n", _my_mpi_id, this->_U, this->_U_hydr, this->_U_stack);

				// we should take care of the external potential here...
				if (this->_external_forces) {
					fact *= exp ((b1 - b2) * (_U_ext - _exchange_energy.U_ext));
				}

				//printf ("(from %d) fact: %lf\n", _my_mpi_id, fact);
				if (drand48() < fact) {
					this->_pt_exchange_accepted ++;
					// send my conf over there
					// store the other guy's
					//printf ("(from %d) accepting exchange...\n", _my_mpi_id);
					PT_serialized_particle_info<number> * buffer_conf;
					PT_energy_info<number> buffer_energy;
					buffer_conf = new PT_serialized_particle_info<number>[this->_N];
					for (int i = 0; i < this->_N; i ++) {
						buffer_conf[i].pos = _exchange_conf[i].pos;
						buffer_conf[i].orientation = _exchange_conf[i].orientation;
					}
					buffer_energy.U = _exchange_energy.U;
					buffer_energy.U_hydr = _exchange_energy.U_hydr;
					buffer_energy.U_stack = _exchange_energy.U_stack;
					buffer_energy.T = _exchange_energy.T;
					buffer_energy.U_ext = _exchange_energy.U_ext;
					buffer_energy.replica_id = _exchange_energy.replica_id;

					_build_exchange_energy ();
					_send_exchange_energy (irresp_id);

					// put my current one in exchange_conf
					_build_exchange_conf ();
					_send_exchange_conf (irresp_id);

					// put the other guy's old conf in my echange vector
					for (int i = 0; i < this->_N; i ++) {
						_exchange_conf[i].pos = buffer_conf[i].pos;
						_exchange_conf[i].orientation = buffer_conf[i].orientation;
					}
					_exchange_energy.U = buffer_energy.U;
					_exchange_energy.U_hydr = buffer_energy.U_hydr;
					_exchange_energy.U_stack = buffer_energy.U_stack;
					_exchange_energy.T = buffer_energy.T;
					_exchange_energy.U_ext = buffer_energy.U_ext;
					_exchange_energy.replica_id = buffer_energy.replica_id;

					// rebuild my conf (which will now become the other guy's
					// old one
					_rebuild_exchange_conf ();
					_rebuild_exchange_energy ();

					delete [] buffer_conf;
				}
				else {
					//printf ("(from %d) rejecting exchange...\n", _my_mpi_id);
					// send back the guy's conf as it was
					_send_exchange_energy (irresp_id);
					_send_exchange_conf (irresp_id);
					// send back either its old one or my old one;
				}
				//fprintf (stderr, "(replica %d) got %lf %lf %lf AAA\n", _my_mpi_id, this->_U, this->_U_hydr, this->_U_stack);
				VMMC_CPUBackend<number>::_compute_energy ();
				//fprintf (stderr, "(replica %d) now %lf %lf %lf AAA\n", _my_mpi_id, this->_U, this->_U_hydr, this->_U_stack);
			}
		}
		else {
			// not responsible. I basically don't know anything, but perhaps
			// more importantly I have no interest in doing so.
			if (resp_id >= 0 && resp_id < _mpi_nprocs) {
				//fprintf (stderr, "(replica %d) had %lf %lf %lf AAA\n", _my_mpi_id, this->_U, this->_U_hydr, this->_U_stack);
				// send my energy and conf
				_build_exchange_energy ();
				_send_exchange_energy (resp_id);

				_build_exchange_conf ();
				_send_exchange_conf (resp_id);

				// wait for the conf to use next;
				_get_exchange_energy (resp_id);
				_get_exchange_conf (resp_id);

				_rebuild_exchange_energy ();
				_rebuild_exchange_conf ();
				//fprintf (stderr, "(replica %d) got %lf %lf %lf AAA\n", _my_mpi_id, this->_U, this->_U_hydr, this->_U_stack);
				//VMMC_CPUBackend<number>::_compute_energy ();
				//fprintf (stderr, "(replica %d) now %lf %lf %lf AAA\n", _my_mpi_id, this->_U, this->_U_hydr, this->_U_stack);
			}
		}
		// debug
		if (fabs(this->_T - _pttemps[_my_mpi_id]) > 1.e-7) {
			fprintf (stderr, "DISASTRO\n\n\n");
		}
		// we should se the forces again if we have swapped conf
		if (this->_external_forces) {
			BaseParticle<number> *p;
			for (int i = 0; i < this->_N; i ++) {
				p = this->_particles[i];
				p->set_ext_potential(curr_step, this->_box_side);
			}
		}

	}

	return;
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_send_exchange_energy (int other_id) {
	//fprintf (stderr, "(replica %d) simian %lf %lf %lf\n", _my_mpi_id, _exchange_energy.U, _exchange_energy.U_hydr, _exchange_energy.U_stack);
	//printf ("(from %d) sending energy info to %d\n", _my_mpi_id, other_id);
	_MPI_send_block_data ((void *)(&_exchange_energy), sizeof(PT_energy_info<number>), other_id);
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_get_exchange_energy (int other_id) {
	//printf ("(from %d) waiting energy info from %d\n", _my_mpi_id, other_id);
	_MPI_receive_block_data ((void *)(&_exchange_energy), sizeof(PT_energy_info<number>), other_id);
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_send_exchange_conf (int other_id) {
	//printf ("(from %d) sending conf info to %d\n", _my_mpi_id, other_id);
	_MPI_send_block_data ((void *)(_exchange_conf), this->_N * sizeof(PT_serialized_particle_info<number>), other_id);
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_get_exchange_conf (int other_id) {
	//printf ("(from %d) waiting conf info from %d\n", _my_mpi_id, other_id);
	_MPI_receive_block_data ((void *)(_exchange_conf), this->_N * sizeof(PT_serialized_particle_info<number>), other_id);
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_build_exchange_conf() {
	for (int i = 0; i < this->_N; i ++) {
		_exchange_conf[i].read_from(this->_particles[i]);
	}
	return;
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_build_exchange_energy() {
	_exchange_energy.U = this->_U;
	_exchange_energy.U_hydr = this->_U_hydr;
	_exchange_energy.U_stack = this->_U_stack;
	_exchange_energy.T = this->_T;
	_exchange_energy.U_ext = _U_ext;
	_exchange_energy.replica_id = _which_replica;
	if (this->_have_us) _exchange_energy.w = this->_w.get_weight (this->_op.get_hb_states(), &(_exchange_energy.weight_index));
}


template <typename number>
void PT_VMMC_CPUBackend<number>::_rebuild_exchange_energy() {
	/*
	this->_U = _exchange_energy.U;
	this->_U_hydr = _exchange_energy.U_hydr;
	this->_U_stack = _exchange_energy.U_stack;
	*/
	//et = this->_U_stack * STCK_FACT_EPS / (STCK_BASE_EPS + STCK_FACT_EPS / b1);
	//fprintf (stderr, "(replica %d) monkey %lf %lf %lf\n", _my_mpi_id, _exchange_energy.U, _exchange_energy.U_hydr, _exchange_energy.U_stack);

	if (_oxDNA2_stacking) this->_U_stack = _exchange_energy.U_stack * (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * this->_T) / (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * _exchange_energy.T);
	else this->_U_stack = _exchange_energy.U_stack * (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * this->_T) / (STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * _exchange_energy.T);
	this->_U = _exchange_energy.U + this->_U_stack - _exchange_energy.U_stack;
	this->_U_hydr = _exchange_energy.U_hydr;
	_which_replica = _exchange_energy.replica_id;
	// N.B. we ignore the temperature, as we need it to decide wether to accept
	// or not, but we don't want to overwrite our current one
}

template <typename number>
void PT_VMMC_CPUBackend<number>::_rebuild_exchange_conf() {
	BaseParticle<number> * p;
	for (int i = 0; i < this->_N; i ++) {
		p = this->_particles[i];
		_exchange_conf[i].write_to(p);
		p->orientationT = p->orientation.get_transpose();
		p->set_positions ();
	}
	VMMC_CPUBackend<number>::_delete_cells ();
	VMMC_CPUBackend<number>::_init_cells ();
	//VMMC_CPUBackend<number>::_update_metainfo ();
	return;
}

template <typename number>
int PT_VMMC_CPUBackend<number>::_MPI_send_block_data (void * data, size_t size, int node_to, int TAG) {
	int ret = MPI_Send((void *)data, size, MPI_CHAR, node_to, TAG, MPI_COMM_WORLD);
	if(ret != MPI_SUCCESS) {
		throw oxDNAException("Error while sending MPI message");
		return -1;
	}
	else {
		return 1;
	}
}

template <typename number>
int PT_VMMC_CPUBackend<number>::_MPI_receive_block_data (void *data, size_t size, int node_from, int TAG) {
	MPI_Status stat;
	int ret = MPI_Recv((void *)data, size, MPI_CHAR, node_from, TAG, MPI_COMM_WORLD, &stat);

	if(ret != MPI_SUCCESS) {
		throw oxDNAException("Error while receving MPI message");
		return -1;
	}
	else {
		return 1;
	}
}

template<typename number>
number PT_VMMC_CPUBackend<number>::get_pt_acc () {
	if (_my_mpi_id == (_mpi_nprocs - 1)) {
		// I am never responsible, so by definition...
		return 0.;
	}
	else {
		return _pt_exchange_accepted / (number) _pt_exchange_tries;
	}
}

template<typename number>
char * PT_VMMC_CPUBackend<number>::get_replica_info_str() {
	sprintf (_replica_info, "%d %d", _my_mpi_id, _which_replica);
	return _replica_info;
}

template class PT_VMMC_CPUBackend<float> ;
template class PT_VMMC_CPUBackend<double> ;


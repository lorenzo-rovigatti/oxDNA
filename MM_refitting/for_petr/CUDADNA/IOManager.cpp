/*
 * IOManager.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include "IOManager.h"
#include "CMMC_CPUBackend.h"
#include "FL_CMMC_CPUBackend.h"
#include "FFS_MD_CPUBackend.h"
#include "VMMC_CPUBackend.h"
#include "PT_VMMC_CPUBackend.h"

IOManager::IOManager() : _debug(false), _restart_step_counter(false) {
	_energy_stream = NULL;
	_log_open = false;
	_log_stream = stderr;
    _allow_log = true;

	sprintf(_log_level_strings[0], "%s", "INFO");
	sprintf(_log_level_strings[1], "%s", "WARNING");
	sprintf(_log_level_strings[2], "%s", "ERROR");
	sprintf(_log_level_strings[3], "%s", "CRITICAL");
	sprintf(_log_level_strings[4], "%s", "DEBUG");
}

IOManager::~IOManager() {
	if(_log_open) fclose(_log_stream);
	if(_energy_stream != NULL) fclose(_energy_stream);
}

// Copy function 

void IOManager::_log(int log_level, const char *format, va_list &ap) {

	if(!_allow_log)
		return;

	if(log_level != LOG_NOTHING) fprintf(_log_stream, "%s: ", _log_level_strings[log_level]);
	vfprintf(_log_stream, format, ap);
	va_end(ap);

	fprintf(_log_stream, "\n");
	fflush(_log_stream);
}

void IOManager::_set_stream(const char *filename) {
	FILE *buff = fopen(filename, "w");
	if(buff == NULL) die("Log file '%s' is not writable", filename);

	_log_stream = buff;
	_log_open = true;
}

void IOManager::log(int log_level, const char *format, ...) {
	va_list ap;
	va_start(ap, format);
	_log(log_level, format, ap);
}

void IOManager::debug(const char *format, ...) {
	if(_debug == true) {
		va_list ap;
		va_start(ap, format);
		_log(LOG_DEBUG, format, ap);
	}
}

void IOManager::die(const char *format, ...) {
	va_list ap;
	va_start(ap, format);
	_log(LOG_CRITICAL, format, ap);

#ifdef HAVE_MPI
	MPI_Abort(MPI_COMM_WORLD,1);
#endif

	fprintf(_log_stream, "#######################################\n######         GAME OVER         ######\n#######################################\n");
	exit(-1);
}

void IOManager::print_version() {
	fprintf(stdout, "oxDNA %d.%d.%d by Lorenzo Rovigatti, Flavio Romano and Petr Sulc (c) 2012\n", VERSION_MAJOR, VERSION_MINOR, VERSION_STAGE);
	exit(-1);
}

void IOManager::get_settings(input_file &inp) {

	char filename[256];
	int tmp;

	if(getInputString(&inp, "log_file", filename, 0) == KEY_FOUND) _set_stream(filename);

	getInputString(&inp, "trajectory_file", _trajectory_file, 1);
	getInputString(&inp, "energy_file", _energy_file, 1);
	if(!(getInputString(&inp, "lastconf_file", _last_conf_file, 0) == KEY_FOUND)) 
                                  strcpy(_last_conf_file,"last_conf.dat");
	
	if(getInputInt(&inp, "print_reduced_conf_every", &tmp, 0) == KEY_FOUND && tmp > 0)
		getInputString(&inp, "reduced_conf_output_dir", _reduced_conf_output_dir, 1);
	if(getInputInt(&inp, "restart_step_counter", &tmp, 0) == KEY_FOUND) _restart_step_counter = (tmp);
	if(getInputInt(&inp, "debug", &tmp, 0) == KEY_FOUND) _debug = (tmp);
	if(getInputInt(&inp, "no_stdout_energy", &tmp, 0) == KEY_FOUND) _no_stdout_energy = (tmp != 0);
	else _no_stdout_energy = false;

#ifdef HAVE_MPI
	char sim_type[256];
	getInputString(&inp, "sim_type", sim_type, 1);
	if (strcmp (sim_type, "PT_VMMC") == 0) {
		char extra[128];
		int myid, nprocs;
		MPI_Comm_rank (MPI_COMM_WORLD, &(myid));
		MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
		sprintf (extra, "%i", myid);
		if (nprocs > 1) {
			strcat (_trajectory_file, extra);
			strcat (_last_conf_file, extra);
			strcat (_energy_file, extra);
		}
	}
#endif
}

void IOManager::init(bool append) {
	ofstream traj_file;

	//if(_restart_step_counter) {
	if (append) {
		_energy_stream = fopen(_energy_file, "w");
		traj_file.open(_trajectory_file);
	}
	else {
		_energy_stream = fopen(_energy_file, "a");
		traj_file.open(_trajectory_file, ios_base::app);
	}

	if(!traj_file.good()) this->die("Trajectory file '%s' is not writable", _trajectory_file);
	traj_file.close();
}

template<typename number>
void IOManager::_print_conf(SimBackend<number> &backend, llint curr_step, ofstream &conf_output, bool reduced) {
	int N = backend._N;

	conf_output.precision(15);

	conf_output << "t = " << curr_step << endl;
	conf_output << "b = " << backend._box_side << " " << backend._box_side << " " << backend._box_side << endl;
	conf_output << "E = " << backend._K/backend._N + backend._U/backend._N << " " << backend._U/backend._N << " " << backend._K/backend._N  << endl;

	for(int i = 0; i < N; i++) {
		LR_matrix<number> oT = backend._particles[i].orientation.get_transpose();
		conf_output << backend._particles[i].pos.x << " " << backend._particles[i].pos.y << " " << backend._particles[i].pos.z << " ";
		conf_output << oT.v1.x << " " << oT.v1.y << " " << oT.v1.z << " ";
		conf_output << oT.v3.x << " " << oT.v3.y << " " << oT.v3.z << " ";
		conf_output << backend._particles[i].vel.x << " " << backend._particles[i].vel.y << " " << backend._particles[i].vel.z << " ";
		conf_output << backend._particles[i].L.x << " " << backend._particles[i].L.y << " " << backend._particles[i].L.z << endl;
	}
}

template<typename number>
void IOManager::print_conf(SimBackend<number> &backend, llint curr_step, bool only_last) {
	//char conf_name[256] = "last_conf.dat";

	ofstream conf_output(_last_conf_file);
	if(!conf_output.good()) this->die("'%s' is not writable", _last_conf_file);
	_print_conf(backend, curr_step, conf_output);
	conf_output.close();

	if(only_last == false) {
		conf_output.open(_trajectory_file, ios_base::app);
		_print_conf(backend, curr_step, conf_output);
		conf_output.close();
	}
}

template<typename number>
void IOManager::print_reduced_conf(SimBackend<number> &backend, llint curr_step) {
	char conf_name[256];

	sprintf(conf_name, "%s/reduced_conf%lld.dat", _reduced_conf_output_dir, curr_step);

	ofstream conf_output(conf_name);
	if(!conf_output.good()) this->die("'%s' is not writable", conf_name);
	//_print_conf(backend, curr_step, conf_output, true);

	conf_output << backend._N_strands << " " << backend._box_side << " " << backend._box_side << " " << backend._box_side << " " << curr_step << endl;

	LR_vector<number> cm(0, 0, 0);
	int counter = 0;
	for(int i = 0; i < backend._N; i++) {
		Particle<number> *p = &backend._particles[i];

		if(p->n3 == P_VIRTUAL) {
			cm = p->pos;
			counter = 1;
		}
		else {
			cm += p->pos;
			counter++;
		}

		if(p->n5 == P_VIRTUAL) {
			cm /= counter;
			conf_output << cm.x << " " << cm.y << " " << cm.z << endl;
		}
	}

	conf_output.close();
}

// print_energy method for MD
template<typename number>
void IOManager::print_energy(MDBackend<number> &backend, llint curr_step) {
		if(_no_stdout_energy == false)
			printf("%15lld % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf\n", curr_step, curr_step*backend._dt, backend._U/backend._N, backend._K/backend._N, backend._K/backend._N + backend._U/backend._N, backend._U_hydr);

		fprintf(_energy_stream, "%lf %lf %lf %lf %lf\n", curr_step*backend._dt, backend._U/backend._N, backend._K/backend._N, backend._K/backend._N + backend._U/backend._N, backend._U_hydr);
		fflush(_energy_stream);
}

// print_energy method for MC
template<typename number>
void IOManager::print_energy(MCBackend<number> &backend, llint curr_step) {
		if(_no_stdout_energy == false) {
			printf("%15lld %10.6lf %10.6lf", curr_step, backend._U/backend._N, backend._U_hydr);
			for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) printf(" %10.6lf", backend._accepted[i]/(float)backend._tries[i]);
			printf("\n");
		}

		fprintf(_energy_stream, "%lld %lf %lf", curr_step, backend._U/backend._N, backend._U_hydr);
		for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) fprintf(_energy_stream, " %lf", backend._accepted[i]/(float)backend._tries[i]);
		fprintf(_energy_stream, "\n");
		fflush(_energy_stream);
}

// print_energy method for CMMC
template<typename number>
void IOManager::print_energy(CMMC_CPUBackend<number> &backend, llint curr_step) {
	if(_no_stdout_energy == false) {
		printf("%15lld %10.6lf %10.6lf", curr_step, backend._U/backend._N, backend._U_hydr);
		for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) printf(" %10.6lf", backend._accepted[i]/(float)backend._tries[i]);

		printf(" %s", backend.get_op_state_str());
		printf ("\n");
	}

	fprintf(_energy_stream, "%lld %lf %lf", curr_step, backend._U/backend._N, backend._U_hydr);
	for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) fprintf(_energy_stream, " %lf", backend._accepted[i]/(float)backend._tries[i]);
	fprintf(_energy_stream, " %s", backend.get_op_state_str());
	fprintf(_energy_stream, "\n");
	fflush(_energy_stream);
}

// print_energy method for VMMC
template<typename number>
void IOManager::print_energy(VMMC_CPUBackend<number> &backend, llint curr_step) {
	if(_no_stdout_energy == false) {
		printf("%15lld %10.6lf %10.6lf", curr_step, backend._U/backend._N, backend._U_hydr);
		for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) printf(" %10.6lf", backend._accepted[i]/(float)backend._tries[i]);

		printf(" %s", backend.get_op_state_str());
		printf ("\n");
	}

	fprintf(_energy_stream, "%lld %lf %lf", curr_step, backend._U/backend._N, backend._U_hydr);
	for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) fprintf(_energy_stream, " %lf", backend._accepted[i]/(float)backend._tries[i]);
	fprintf(_energy_stream, " %s", backend.get_op_state_str());
	fprintf(_energy_stream, "\n");
	fflush(_energy_stream);
}

// print_energy method for PT_VMMC
#ifdef HAVE_MPI
template<typename number>
void IOManager::print_energy(PT_VMMC_CPUBackend<number> &backend, llint curr_step) {
	if(_no_stdout_energy == false) {
		printf("%15lld %10.6lf %10.6lf", curr_step, backend._U/backend._N, backend._U_hydr);
		for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) printf(" %10.6lf", backend._accepted[i]/(float)backend._tries[i]);

		printf (" %5.3lf", backend.get_pt_acc ());

		printf(" %s", backend.get_op_state_str());
		printf (" REPLICA %d", backend.get_mpi_id ());
		printf (" CONF %d", backend.get_which_replica ());
		printf ("\n");
	}

	fprintf(_energy_stream, "%lld %lf %lf", curr_step, backend._U/backend._N, backend._U_hydr);
	for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) fprintf(_energy_stream, " %lf", backend._accepted[i]/(float)backend._tries[i]);
	fprintf(_energy_stream, " %s", backend.get_op_state_str());
	fprintf(_energy_stream, " %s", backend.get_replica_info_str());
	fprintf(_energy_stream, "\n");
	fflush(_energy_stream);
}
#endif

// print_energy method for CMMC
template<typename number>
void IOManager::print_energy(FL_CMMC_CPUBackend<number> &backend, llint curr_step) {
	if(_no_stdout_energy == false) {
		printf("%15lld %10.6lf %10.6lf", curr_step, backend._U/backend._N, backend._U_hydr);
		for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) printf(" %10.6lf", backend._accepted[i]/(float)backend._tries[i]);

		printf(" %s", backend.get_op_state_str());
		printf(" %s", backend.get_fl_str());
		printf(" %lf", backend.get_e_stack ());
		printf ("\n");
	}

	fprintf(_energy_stream, "%lld %lf %lf", curr_step, backend._U/backend._N, backend._U_hydr);
	for(int i = 0; i < MC_MOVES; i++) if(backend._tries[i] > 0) fprintf(_energy_stream, " %lf", backend._accepted[i]/(float)backend._tries[i]);
	fprintf(_energy_stream, " %s", backend.get_op_state_str());
	fprintf(_energy_stream, " %s", backend.get_fl_str());
	fprintf(_energy_stream, " %lf", backend.get_e_stack ());
	fprintf(_energy_stream, "\n");
	fflush(_energy_stream);
}

template<typename number>
void IOManager::print_energy(FFS_MD_CPUBackend<number> &backend, llint curr_step) {
	if(_no_stdout_energy == false) {
		printf("%15lld % 10.6lf % 10.6lf % 10.6lf % 10.6lf % 10.6lf ", curr_step, curr_step*backend._dt, backend._U/backend._N, backend._K/backend._N, backend._K/backend._N + backend._U/backend._N, backend._U_hydr);
		printf(" %s", backend.get_op_state_str());
		printf ("\n");
	}

	fprintf(_energy_stream, "%lf %lf %lf %lf %lf ", curr_step*backend._dt, backend._U/backend._N, backend._K/backend._N, backend._K/backend._N + backend._U/backend._N, backend._U_hydr);
	fprintf(_energy_stream, " %s", backend.get_op_state_str());
	fprintf(_energy_stream, "\n");
	fflush(_energy_stream);
}
template void IOManager::_print_conf(SimBackend<float> &backend, llint curr_step, ofstream &, bool);
template void IOManager::_print_conf(SimBackend<double> &backend, llint curr_step, ofstream &, bool);

template void IOManager::print_conf(SimBackend<float> &, llint, bool);
template void IOManager::print_conf(SimBackend<double> &, llint, bool);

template void IOManager::print_reduced_conf(SimBackend<float> &, llint);
template void IOManager::print_reduced_conf(SimBackend<double> &, llint);

template void IOManager::print_energy(MDBackend<float> &, llint);
template void IOManager::print_energy(MDBackend<double> &, llint);

template void IOManager::print_energy(MCBackend<float> &, llint);
template void IOManager::print_energy(MCBackend<double> &, llint);

template void IOManager::print_energy(FFS_MD_CPUBackend<double> &, llint);
template void IOManager::print_energy(FFS_MD_CPUBackend<float> &, llint);

template void IOManager::print_energy(CMMC_CPUBackend<double> &, llint);
template void IOManager::print_energy(CMMC_CPUBackend<float> &, llint);

template void IOManager::print_energy(VMMC_CPUBackend<double> &, llint);
template void IOManager::print_energy(VMMC_CPUBackend<float> &, llint);

#ifdef HAVE_MPI
template void IOManager::print_energy(PT_VMMC_CPUBackend<double> &, llint);
template void IOManager::print_energy(PT_VMMC_CPUBackend<float> &, llint);
#endif

template void IOManager::print_energy(FL_CMMC_CPUBackend<double> &, llint);
template void IOManager::print_energy(FL_CMMC_CPUBackend<float> &, llint);

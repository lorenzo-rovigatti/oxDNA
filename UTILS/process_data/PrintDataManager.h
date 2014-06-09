/*
 * PrintDataManager.h
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#ifndef PRINTMANAGER_H_
#define PRINTMANAGER_H_

#include <cstring>
#include <ctime>

extern "C" {
#include "../../timing/timing.h"
#include "../../time_scales/time_scales.h"
}

#include "../../defs.h"

#include "../../IOManager.h"
#include "ProcessData_Backend.h"

struct double4;
struct float4;

class PrintDataManager {
protected:
	IOManager &_IO;
	ProcessData_Backend<double> *_backend;
	input_file _input;
	time_scale _time_scale_manager;
	int _time_scale;
	llint _steps;
	llint _cur_step;
	llint _start_step;
	llint _max_steps;
	int _seed;
	int _print_reduced_conf_every;
	int _print_energy_every;
	int _restart_step_counter;
	int _iteration;
	char _conf_file[256];
	void _get_options();

	ifstream conf_input;

public:
	PrintDataManager(IOManager &IO, const char *input_file,int confs_to_skip=0);
	virtual ~PrintDataManager();
	
	static bool stop;
	void load_options();
	void init(const char* input_fname);
	void run();
	void clean();
	ostream& OutputState(ostream& out,int t=0); //outputs energies between particles

};

void gbl_termhandler (int);

#endif /* SIMMANAGER_H_ */

/**
 * @file    MC_CPUBackend2.h
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#ifndef MC_CPUBACKEND2_H_
#define MC_CPUBACKEND2_H_

#include "MCBackend.h"
#include "MCMoves/BaseMove.h"
#include "MCMoves/MoveFactory.h"

/**
 * @brief Manages a MC simulation on CPU. It supports NVT and NPT simulations
 */

class MC_CPUBackend2: public MCBackend {
protected:
	std::vector<MovePtr> _moves;
	int _N_moves;
	std::shared_ptr<ConfigInfo> _config_info;
	std::string _info_str;
	void _compute_energy(); // to silence compiler, for now
	number _accumulated_prob;

public:
	MC_CPUBackend2();
	virtual ~MC_CPUBackend2();

	virtual void get_settings(input_file &inp);
	void init();

	void sim_step();
	void add_move(std::string move_string, input_file &sim_inp);

	void print_observables();

	virtual void print_equilibration_info();
};

#endif /* MC_CPUBACKEND2_H_ */

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
template<typename number>
class MC_CPUBackend2: public MCBackend<number> {
protected:
	std::vector <BaseMove<number> *> _moves;
	int _N_moves;
	ConfigInfo<number> _MC_Info;
	std::string _info_str;
	void _compute_energy(); // to silence compiler, for now
	number _accumulated_prob;
	number _verlet_skin;

public:
	MC_CPUBackend2();
	virtual ~MC_CPUBackend2();

	virtual void get_settings(input_file &inp);
	void init();

	void sim_step(llint cur_step);
	void add_move (std::string move_string, input_file &sim_inp);
	
	void print_observables(llint curr_step); 
};

#endif /* MC_CPUBACKEND2_H_ */

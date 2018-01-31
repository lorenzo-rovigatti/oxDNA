/**
 * @file    MD_CPUBackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MD_CPUBACKEND_H_
#define MD_CPUBACKEND_H_

#include "MDBackend.h"
#include "MCMoves/VolumeMove.h"

template <typename number> class BaseThermostat;

/**
 * @brief Manages a MD simulation on CPU. It supports NVE and NVT simulations
 */
template<typename number>
class MD_CPUBackend: public MDBackend<number> {
protected:
	BaseThermostat<number> * _thermostat;
	BaseMove<number> *_V_move;
	bool _compute_stress_tensor;
	int _stress_tensor_avg_every;
	int _stress_tensor_counter;
	LR_matrix<double> _stress_tensor;

	void _first_step(llint cur_step);
	void _compute_forces();
	void _second_step();

	void _update_forces_and_stress_tensor(BaseParticle<number> *p, BaseParticle<number> *q);
	void _update_kinetic_stress_tensor(BaseParticle<number> *p);
	void _update_backend_info();

public:
	MD_CPUBackend();
	virtual ~MD_CPUBackend();

	void init();
	void get_settings (input_file &inp);
	void sim_step(llint cur_step);
	void activate_thermostat();
};

#endif /* MD_CPUBACKEND_H_ */

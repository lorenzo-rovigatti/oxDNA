/**
 * @file    MCBackend.h
 * @date    25/nov/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MCBACKEND_H_
#define MCBACKEND_H_

#define MC_ENSEMBLE_NVT 0
#define MC_ENSEMBLE_NPT 1

#define MC_MOVES 4
#define MC_MOVE_TRANSLATION 0
#define MC_MOVE_ROTATION 1
#define MC_MOVE_VOLUME 2
#define MC_MOVE_CLUSTER_SIZE 3

#include "SimBackend.h"

/**
 * @brief Abstract class for MC simulations.
 *
 * This class sets up a basic MC simulation. It does not do any sensible computation but
 * takes care of the most basic input/output operations associated with MC simulations.
 *
 * @verbatim
 ensemble = nvt|npt (ensemble of the simulation)
 [check_energy_every = <float> (oxDNA will compute the energy from scratch, compare it with the current energy and throw an error if the difference is larger then check_energy_threshold. Defaults to 10.)]
 [check_energy_threshold = <float> (threshold for the energy check. Defaults to 1e-2f for single precision and 1e-6 for double precision.)]

 delta_translation = <float> (controls the trial translational displacement, which is a randomly chosen number between -0.5*delta and 0.5*delta for each direction.)
 delta_rotation = <float> (controls the angular rotational displacement, given by a randomly chosen angle between -0.5*delta and 0.5*delta radians.)
 delta_volume = <float> (controls the volume change in npt simulations.)

 P = <float> (the pressure of the simulation. Used only if ensemble == npt.)

 [adjust_moves = <bool> (if true, oxDNA will run for equilibration_steps time steps while changing the delta of the moves in order to have an optimal acceptance ratio. It does not make sense if equilibration_steps is 0 or not given. Defaults to false)]
 @endverbatim
 */

class MCBackend: public SimBackend {
protected:
	int _MC_moves;
	bool _overlap;
	bool _adjust_moves;
	std::vector<number> _delta;
	std::vector<llint> _tries;
	std::vector<llint> _accepted;
	llint _MC_equilibration_steps;
	int _ensemble;
	int _check_energy_every;
	int _check_energy_counter;
	number _check_energy_threshold;

	virtual void _compute_energy() = 0;

public:
	MCBackend();
	virtual ~MCBackend();

	void get_settings(input_file &inp);

	virtual void print_observables();
};

#endif /* MCBACKEND_H_ */


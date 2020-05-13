/*
 * Bussi thermostat class.
 *
 *  Created on: May 07, 2015
 *      Author: Flavio
 */

#ifndef BUSSI_THERMOSTAT_
#define BUSSI_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief This class implements the thermostat introduced in the paper
 * "Canonical sampling through velocity rescaling"
 * by Bussi, Donadio and Parrinello (JCP 126, 014101, 2007).
 * This is a canonical thermostat which is good for atomic and molecular systems. It contains
 * a stochastic evolution of the kinetic energy. It is supposed to be an improvement over 
 * Nose-Hoover or other thermostats. In particular, it is able to sample the fluctuations 
 * of the kinetic energy. Do NOT use for DNA or other coarse grained models since it does
 * not take the solvent into account.
 *
 * @verbatim
 newtonian_steps = <int> (number of integration timesteps after which the thermostat acts. Can be 1.)
 bussi_tau = <int> (correlation time, in time steps, for the stochastic evolution of the kinetic energy)
 @endverbatim
 */

class BussiThermostat: public BaseThermostat {
protected:
	/// interval, in time steps, over which the dynamics is ballistic.
	int _newtonian_steps;

	/// correlation time for the stochastic dynamics of the kinetic energy (in time steps)
	int _tau;

	/// number used in the Bussi formula to integrate the kinetic energy
	number _exp_dt_tau;

	/// stored value for the stochastic dynamic of the kinetic energy associated to translations
	number _K_t;

	/// stored value for the stochastic dynamic of the kinetic energy associated to rotations
	number _K_r;

	void _update_K(number &K, int degrees_of_freedom);

	// Bussi's methods
	number _sum_noises(int nn);
	number _gamdev(const int ia);
	int _current_translational_degrees_of_freedom();
	int _current_rotational_degrees_of_freedom();

public:
	BussiThermostat();
	virtual ~BussiThermostat();

	void get_settings(input_file &inp);
	void init();
	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif // REFRESH_THERMOSTAT_

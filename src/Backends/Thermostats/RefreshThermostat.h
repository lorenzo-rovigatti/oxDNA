/*
 * Refresh thermostat class.
 */

#ifndef REFRESH_THERMOSTAT_
#define REFRESH_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief Incapsulates a simple Brownian thermostat that refreshes each particle's momenta every n timesteps.
 *
 * @verbatim
 newtonian_steps = <int> (number of integration timesteps after which momenta are refreshed)
 @endverbatim
 */

class RefreshThermostat: public BaseThermostat {
private:
	int _newtonian_steps;
	number _rescale_factor;
public:
	RefreshThermostat();
	virtual ~RefreshThermostat();

	void get_settings(input_file &inp);
	void init();
	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif // REFRESH_THERMOSTAT_

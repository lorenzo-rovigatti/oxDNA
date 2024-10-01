/*
 * No (fake) thermostat class.
 */

#ifndef NO_THERMOSTAT_
#define NO_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief Incapsulates a lack of thermostat, if that makes any sense.
 */

class NoThermostat : public BaseThermostat {
public:
	NoThermostat();
	virtual ~NoThermostat();

	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif // NO_THERMOSTAT_

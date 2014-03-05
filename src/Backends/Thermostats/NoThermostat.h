/*
 * No (fake) thermostat class.
 */

#ifndef NO_THERMOSTAT_
#define NO_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief Incapsulates a lack of thermostat, if that makes any sense.
 */
template<typename number>
class NoThermostat : public BaseThermostat<number> {
public:
	NoThermostat();
	virtual ~NoThermostat();

	void apply(BaseParticle<number> **particles, llint curr_step);
};

#endif // NO_THERMOSTAT_

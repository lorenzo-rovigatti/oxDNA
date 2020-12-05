#ifndef BASE_THERMOSTAT_
#define BASE_THERMOSTAT_

#include "../../Utilities/Utils.h"
#include "../../Utilities/oxDNAException.h"
#include "../../defs.h"
#include "../../Particles/BaseParticle.h"

/**
 * @brief This class is the basic thermostat interface.
 */

class IBaseThermostat {
public:
	virtual ~IBaseThermostat() {
	}
	;

	/**
	 * @brief function that gets the settings from the input file.
	 * Sets the temperature from the input file.
	 *
	 * @param inp input file object reference
	 */
	virtual void get_settings(input_file &inp) = 0;

	/**
	 * @brief init function for the thermostats.
	 */
	virtual void init() = 0;
};

/**
 * @brief This class is meant to be used so that all other thermostat classes can inherit
 * from here
 */

class BaseThermostat: public virtual IBaseThermostat {
protected:
	number _T;
	bool _supports_shear;
	bool _lees_edwards;
	number _shear_rate;

	virtual void _on_T_update();

public:
	BaseThermostat();
	virtual ~BaseThermostat() {
	}

	virtual void get_settings(input_file &inp);
	virtual void init() {

	}

	/**
	 * @brief this method is what the MD_CPUBackend calls to apply the
	 * thermostat to the system
	 *
	 * This method is purely virtual, meaning that it has to be implemented
	 * in all the classes that inherit it.
	 *
	 * @param particles array of BaseParticle objects
	 * @param curr_step current step of the simulation.
	 */
	virtual void apply(std::vector<BaseParticle *> &particles, llint curr_step) = 0;
};

using ThermostatPtr = std::shared_ptr<BaseThermostat>;

#endif // BASE_THERMOSTAT_

#ifndef BASE_THERMOSTAT_
#define BASE_THERMOSTAT_

#include "../../Utilities/Utils.h"
#include "../../Utilities/oxDNAException.h"
#include "../../defs.h"
#include "../../Particles/BaseParticle.h"

/**
 * @brief This class is the basic thermostat interface.
 */
template <typename number>
class IBaseThermostat {
public:
	virtual ~IBaseThermostat() {};

	/**
	 * @brief function that gets the settings from the input file.
	 * Sets the temperature from the input file.
	 *
	 * @param inp input file object reference
	 */
	virtual void get_settings(input_file &inp) = 0;

	/**
	 * @brief init function for the thermostats.
	 *
	 * @param N number of particles
	 */
	virtual void init(int N) = 0;
};

/**
 * @brief This class is meant to be used so that all other thermostat classes can inherit
 * from here
 */
template <typename number>
class BaseThermostat : public virtual IBaseThermostat<number> {
protected:
	int _N_part;
	number _T;

public:
	BaseThermostat () : _N_part(0), _T ((number) 0.f){}
	virtual ~BaseThermostat () {}

	virtual void get_settings (input_file &inp);
	virtual void init(int N) { _N_part = N; }

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
	virtual void apply(BaseParticle<number> **particles, llint curr_step) = 0;
};

template<typename number>
void BaseThermostat<number>::get_settings(input_file &inp) {
	// we get the temperature from the input file;
	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number> (raw_T);
}

#endif // BASE_THERMOSTAT_


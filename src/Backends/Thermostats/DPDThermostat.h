#ifndef DPD_THERMOSTAT_
#define DPD_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief Implementation of the DPD thermostat
 *
 * the DPD thermostat adds a random force and a relative drag force to
 * the velocities and angular momenta of the particles. It requires dt and
 * the friction coefficient to be specified in the input file.
 * See Peters, Europhys. Lett. 66 (3), pp. 311-317 (2004).
 *
 * @verbatim
 DPD_zeta = <float> (translational damping coefficient for the DPD thermostat.)
 @endverbatim
 */

class DPDThermostat: public BaseThermostat {
protected:
	/// Integration time step
	number _dt, _sqrt_dt;

	/// Translational friction coefficient
	number _zeta;

	/// DPD cut-off
	number _rcut;

	/// Reduced mass of a particle pair
	number _reduced_mass;

	/// Exponent of the dissipative weight function
	number _exponent;

public:
	DPDThermostat();
	virtual ~DPDThermostat();

	void get_settings(input_file &inp);
	void init();
	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif // DPD_THERMOSTAT_

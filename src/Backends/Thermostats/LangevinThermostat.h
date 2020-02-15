#ifndef LANGEVIN_THERMOSTAT_
#define LANGEVIN_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief Implementation of the Langevin thermostat
 *
 * the Langevin thermostat adds random force and a drag force to the velocities
 * and angular momenta of the particles. It requires dt and diffusion
 * coefficient to be specified in the input file.
 *
 * @verbatim
 gamma_trans = <float> (translational damping coefficient for the Langevin thermostat. Either this or diff_coeff should be specified in the input file.)
 @endverbatim
 */

class LangevinThermostat: public BaseThermostat {
protected:
	/// Integration time step
	number _dt;

	/// Translational diffusion coefficient
	number _diff_coeff_trans;

	/// Rotational diffuction coefficient = diff_coeff_trans * 3.
	number _diff_coeff_rot;

	/// Defined as sqrt(2 _diff_coeff_trans  / _dt)
	number _rescale_factor_trans;

	/// Defined as sqrt(2 _diff_coeff_rot  / _dt)
	number _rescale_factor_rot;

	/// Translational damping coeff = diff_coeff_trans / T
	number _gamma_trans;

	/// Angular velocity damping coefficient = diff_coeff_rot / T
	number _gamma_rot;

public:
	LangevinThermostat();
	virtual ~LangevinThermostat();

	void get_settings(input_file &inp);
	void init();
	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif // LANGEVIN_THERMOSTAT_

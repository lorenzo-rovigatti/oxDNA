#ifndef BROWNIAN_THERMOSTAT_
#define BROWNIAN_THERMOSTAT_

#include "BaseThermostat.h"

/**
 * @brief Implementation of the Brownian thermostat
 *
 * the brownian thermostat refreshes the velocities and angular momenta of
 * a fraction of the particles at regular intervals. The fraction of
 * particles that have their velocities (and angular momenta) refreshed
 * determines how strong the thermostat is.
 *
 * @verbatim
newtonian_steps = <int> (number of integration timesteps after which momenta are refreshed)
pt = <float> (probability of refreshing the momenta of each particle)
diff_coeff = <float> (base diffusion coefficient. Either pt or diff_coeff should be specified in the input file)
@endverbatim
 */

class BrownianThermostat : public BaseThermostat {
protected:
	int _newtonian_steps;
	number _pt, _pr, _dt;
	number _diff_coeff;
	number _rescale_factor;
public:
	BrownianThermostat();
	virtual ~BrownianThermostat();

	void get_settings(input_file &inp);
	void init();
	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif // BROWNIAN_THERMOSTAT_

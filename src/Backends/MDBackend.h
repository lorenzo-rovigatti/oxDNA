/**
 * @file    MDBackend.h
 * @date    23/nov/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MDBACKEND_H_
#define MDBACKEND_H_

#include "SimBackend.h"

using namespace std;

/**
 * @brief Abstract class for MD simulations.
 *
 * This class sets up a basic MD simulation. It does not do any sensible computation but
 * takes care of the most basic input/output operations associated with MD simulations.
 *
 * @verbatim
[reset_initial_com_momentum = <bool> (if true the momentum of the centre of mass of the initial configuration will be set to 0. Defaults to false to enforce the reproducibility of the trajectory)]
[reset_com_momentum = <bool> (if true the momentum of the centre of mass will be set to 0 each time fix_diffusion is performed. Defaults to false to enforce the reproducibility of the trajectory)]
@endverbatim
 */
template<typename number>
class MDBackend: public SimBackend<number> {
protected:
	number _dt;
	bool _refresh_velocities;
	bool _reset_initial_com_momentum;
	bool _reset_com_momentum;

	// timers
	Timer *_timer_first_step;
	Timer *_timer_forces;
	Timer *_timer_lists;
	Timer *_timer_thermostat;

	/**
	 * @brief Set the momentum of the centre of mass to 0.
	 */
	virtual void _reset_momentum();

	void _generate_vel();

public:
	MDBackend();
	virtual ~MDBackend();

	void get_settings(input_file &inp);
	void init();
	void fix_diffusion();
};

#endif /* MDBACKEND_H_ */

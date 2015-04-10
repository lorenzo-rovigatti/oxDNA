/**
 * @file    SawtoothForce.h
 * @date    8/apr/2015
 * @author  Megan
 *
 */

#ifndef SAWTOOTHFORCE_H_
#define SAWTOOTHFORCE_H_

#include "BaseForce.h"

/**
 * @brief Force that increases stepwise
 * 
 * This is a simple force, where a particle is pulled with a force that 
 * grows stepwise with time, jumping up by increment every wait_time steps, and which acts in a fixed direction. 
 * The plot of force vs time looks like a sawtooth or a staircase.
 *
 * All units are oxDNA units (time, force = energy / distance, ...)
 *
 * arguments:
@verbatim
particle = <int> (particle to apply the force to. -1 applies it to all particles.)
F0 = <float> (Initial force)
wait_time = <float> (time interval over which the force is constant. Units are (MD/MC) steps.)
increment = <float> (amount by which to increment the force every wait_time steps.)
@endverbatim
 */
template<typename number>
class SawtoothForce : public BaseForce<number> {
private:
	int _particle;
	float _wait_time;
	float _increment;

public:
	SawtoothForce();
	virtual ~SawtoothForce();

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential (llint step, LR_vector<number> &pos);
};

#endif /* SAWTOOTHFORCE_H_ */

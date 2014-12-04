/**
 * @file    ConstantRateForce.h
 * @date    18/oct/2011
 * @author  Flavio
 *
 */

#ifndef CONSTANTRATEFORCE_H_
#define CONSTANTRATEFORCE_H_

#include "BaseForce.h"

/**
 * @brief Force that increases with a constant (possibly 0) rate
 *
 * This is a simple force, where a particle is pulled with a force that 
 * is either constant or grows linearly with time in a fixed direction.
 * All units are oxDNA units (time, force = energy / distance, ...)
 *
 * arguments:
@verbatim
particle <int> particle to apply the force to. -1 applies it to all particles.
F0 <float> Initial force
rate <float> growth rate of the force. It is [oxDNA energy units / (oxDNA distance units * (MD/MC) steps] per 
@endverbatim
 */
template<typename number>
class ConstantRateForce : public BaseForce<number> {
private:
	int _particle;

public:
	ConstantRateForce();
	virtual ~ConstantRateForce();

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential (llint step, LR_vector<number> &pos);
};

#endif /* CONSTANTRATEFORCE_H_ */

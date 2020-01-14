/**
 * @file    LJCone.h
 * @date    24/aug/2017
 * @author  Lorenzo
 *
 */

#ifndef REPULSIVE_CONE_H_
#define REPULSIVE_CONE_H_

#include "BaseForce.h"

/// pos0 * (x,y,z) + position = 0 is the definition of the plane.
/// The pos0 vector is pointing to the halfplane where the
/// repulsion is not acting!
class LJCone: public BaseForce {
private:
	int _particle;
	bool _only_repulsive;
	bool _generate_inside;
	BaseBox *_box;

public:
	int _n;
	number _alpha;
	number _sin_alpha, _cos_alpha, _tan_alpha;
	number _sigma;
	number _cutoff;

	LJCone();
	virtual ~LJCone() {
	}

	void get_settings(input_file &);
	void init(std::vector<BaseParticle *> &, int, BaseBox *);

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // REPULSIVE_CONE_H_

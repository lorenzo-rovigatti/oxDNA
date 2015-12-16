/**
 * @file    LJWall.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef LJWALL_H_
#define LJWALL_H_

#include "BaseForce.h"

/// pos0 * (x,y,z) + position = 0 is the definition of the plane.
/// The pos0 vector is pointing to the halfplane where the
/// repulsion is not acting!
template<typename number>
class LJWall : public BaseForce<number> {
private:
	int _particle;
	bool _only_repulsive;
	bool _generate_inside;

public:
	int _n;
	number _position;
	number _sigma;
	number _cutoff;

	LJWall ();
	virtual ~LJWall() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // LJWALL_H

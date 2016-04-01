/**
 * @file    AlignmentField.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef ALIGNEMENT_FIED_H
#define ALIGNEMENT_FIED_H

#include "BaseForce.h"

/// pos0 * (x,y,z) + position = 0 is the definition of the plane.
/// The pos0 vector is pointing to the halfplane where the
/// repulsion is not acting!
template<typename number>
class AlignmentField : public BaseForce<number> {
private:
	int _particle;
	number _F;
	int _v_idx;
	LR_vector<number> * _v_ptr;

public:
	int _n;
	number _position;
	number _sigma;
	number _cutoff;

	AlignmentField ();
	virtual ~AlignmentField() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // ALIGNEMENT_FIED_H

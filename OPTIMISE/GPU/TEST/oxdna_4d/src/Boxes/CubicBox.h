/*
 * CubicBox.h
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#ifndef CUBICBOX_H_
#define CUBICBOX_H_

#include "BaseBox.h"

/**
 * @brief A cubic simulation box.
 */

class CubicBox: public BaseBox {
protected:
	number _side;
	LR_vector _sides;

public:
	CubicBox();
	virtual ~CubicBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly, number Lz);

	virtual LR_vector min_image(const LR_vector &v1, const LR_vector &v2) const;
	virtual number sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const;

	virtual LR_vector normalised_in_box(const LR_vector &v);
	LR_vector box_sides() const override;
	virtual number V() { return _side*_side*_side; }
	
	virtual LR_vector get_abs_pos(BaseParticle * p); 
	virtual void shift_particle (BaseParticle * p, LR_vector &amount);
};

#endif /* CUBICBOX_H_ */

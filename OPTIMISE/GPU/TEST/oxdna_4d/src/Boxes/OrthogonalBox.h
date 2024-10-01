/*
 * OrthogonalBox.h
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#ifndef ORTHOGONALBOX_H_
#define ORTHOGONALBOX_H_

#include "BaseBox.h"

/**
 * @brief A cubic simulation box.
 */

class OrthogonalBox: public BaseBox {
protected:
	LR_vector _sides;

public:
	OrthogonalBox();
	virtual ~OrthogonalBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly, number Lz);

	virtual LR_vector min_image(const LR_vector &v1, const LR_vector &v2) const;
	virtual number sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const;

	virtual LR_vector normalised_in_box(const LR_vector &v);
	LR_vector box_sides() const override;
	virtual number V() { return _sides.x*_sides.y*_sides.z; }

	virtual LR_vector get_abs_pos(BaseParticle * p); 
	virtual void shift_particle (BaseParticle *p, LR_vector &amount);
};

#endif /* ORTHOGONALBOX_H_ */

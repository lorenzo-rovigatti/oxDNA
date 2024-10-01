/*
 * LeesEdwardsCubicBox.h
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#ifndef SRC_BOXES_LEESEDWARDSCUBICBOX_H_
#define SRC_BOXES_LEESEDWARDSCUBICBOX_H_

#include "CubicBox.h"


class LeesEdwardsCubicBox: public CubicBox {
protected:
	number _factor;

public:
	LeesEdwardsCubicBox();
	virtual ~LeesEdwardsCubicBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly, number Lz);

	virtual LR_vector min_image(const LR_vector &v1, const LR_vector &v2) const;
	virtual number sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const;

//	virtual void shift_particle (BaseParticle *p, LR_vector &amount);
};

#endif /* SRC_BOXES_LEESEDWARDSCUBICBOX_H_ */

/*
 * LeesEdwardsCubicBox.h
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#ifndef SRC_BOXES_LEESEDWARDSCUBICBOX_H_
#define SRC_BOXES_LEESEDWARDSCUBICBOX_H_

#include "CubicBox.h"

template<typename number>
class LeesEdwardsCubicBox: public CubicBox<number> {
protected:
	number _factor;

public:
	LeesEdwardsCubicBox();
	virtual ~LeesEdwardsCubicBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly, number Lz);

	virtual LR_vector<number> min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const;
	virtual number sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const;

//	virtual void shift_particle (BaseParticle<number> *p, LR_vector<number> &amount);
};

#endif /* SRC_BOXES_LEESEDWARDSCUBICBOX_H_ */

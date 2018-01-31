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
template<typename number>
class CubicBox: public BaseBox<number> {
protected:
	number _side;
	LR_vector<number> _sides;

public:
	CubicBox();
	virtual ~CubicBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly, number Lz);

	virtual LR_vector<number> min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const;
	virtual number sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const;

	virtual LR_vector<number> normalised_in_box(const LR_vector<number> &v);
	virtual LR_vector<number> &box_sides() ;
	virtual number V() { return _side*_side*_side; }

	virtual void apply_boundary_conditions(BaseParticle<number> **particles, int N);
	
	virtual LR_vector<number> get_abs_pos(BaseParticle<number> * p); 
	virtual void shift_particle (BaseParticle<number> * p, LR_vector<number> &amount);
};

#endif /* CUBICBOX_H_ */

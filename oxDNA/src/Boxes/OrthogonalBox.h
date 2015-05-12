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
template<typename number>
class OrthogonalBox: public BaseBox<number> {
protected:
	LR_vector<number> _sides;

public:
	OrthogonalBox();
	virtual ~OrthogonalBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly, number Lz);

	virtual LR_vector<number> min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const;
	virtual number sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const;

	virtual LR_vector<number> normalised_in_box(const LR_vector<number> &v);
	virtual LR_vector<number> &box_sides() ;

	virtual void apply_boundary_conditions(BaseParticle<number> **particles, int N);
};

#endif /* ORTHOGONALBOX_H_ */

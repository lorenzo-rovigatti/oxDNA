/*
 * BaseBox.h
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#ifndef BASEBOX_H_
#define BASEBOX_H_

#include "../defs.h"
#include "../Particles/BaseParticle.h"
#include "../Utilities/parse_input/parse_input.h"

/**
 * @brief Abstract class defining a generic simulation box.
 */
template<typename number>
class BaseBox {
public:
	BaseBox();
	virtual ~BaseBox();

	virtual void get_settings(input_file &inp) = 0;
	virtual void init(number Lx, number Ly, number Lz) = 0;

	/**
	 * @brief Returns the minimum image distance between the given two vectors.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual LR_vector<number> min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const = 0;
	/**
	 * @brief Returns the minimum image distance between the two given particles.
	 *
	 * @param p
	 * @param q
	 * @return
	 */
	virtual LR_vector<number> min_image(const BaseParticle<number> &p, const BaseParticle<number> &q) {
		return min_image(p.pos, q.pos);
	}
	/**
	 * @brief Returns the minimum image distance between the two given particles.
	 *
	 * @param p
	 * @param q
	 * @return
	 */
	virtual LR_vector<number> min_image(const BaseParticle<number> *p, const BaseParticle<number> *q) {
		return min_image(p->pos, q->pos);
	}

	/**
	 * @brief Returns the square of the minimum image distance between the two given vectors.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual number sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const = 0;
	/**
	 * @brief Returns the square of the minimum image distance between the two given particles.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual number sqr_min_image_distance(const BaseParticle<number> &p, const BaseParticle<number> &q) {
		return sqr_min_image_distance(p.pos, q.pos);
	}
	/**
	 * @brief Returns the square of the minimum image distance between the two given particles.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual number sqr_min_image_distance(const BaseParticle<number> *p, const BaseParticle<number> *q) {
		return sqr_min_image_distance(p->pos, q->pos);
	}

	/**
	 * @brief Brings back v in the box and returns its normalised components.
	 *
	 * @param v
	 * @return
	 */
	virtual LR_vector<number> normalised_in_box(const LR_vector<number> &v) = 0;

	/**
	 * @brief Returns the length of the box sides.
	 *
	 * @return
	 */
	virtual LR_vector<number> &box_sides() = 0;

	/**
	 * @brief Currently not implemented.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void apply_boundary_conditions(BaseParticle<number> **particles, int N) = 0;
};

#endif /* BASEBOX_H_ */

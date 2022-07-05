/*
 * BaseBox.h
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#ifndef BASEBOX_H_
#define BASEBOX_H_

#include "../defs.h"
#include "../Utilities/parse_input/parse_input.h"

class BaseParticle;

/**
 * @brief Abstract class defining a generic simulation box.
 */
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
	virtual LR_vector min_image(const LR_vector &v1, const LR_vector &v2) const = 0;
	/**
	 * @brief Returns the minimum image distance between the two given particles.
	 *
	 * @param p
	 * @param q
	 * @return
	 */
	virtual LR_vector min_image(const BaseParticle &p, const BaseParticle &q);
	/**
	 * @brief Returns the minimum image distance between the two given particles.
	 *
	 * @param p
	 * @param q
	 * @return
	 */
	virtual LR_vector min_image(const BaseParticle *p, const BaseParticle *q);

	/**
	 * @brief Returns the square of the minimum image distance between the two given vectors.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual number sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const = 0;
	/**
	 * @brief Returns the square of the minimum image distance between the two given particles.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual number sqr_min_image_distance(const BaseParticle &p, const BaseParticle &q);

	/**
	 * @brief Returns the square of the minimum image distance between the two given particles.
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	virtual number sqr_min_image_distance(const BaseParticle *p, const BaseParticle *q);

	/**
	 * @brief Brings back v in the box and returns its normalised components.
	 *
	 * @param v
	 * @return
	 */
	virtual LR_vector normalised_in_box(const LR_vector &v) = 0;

	/**
	 * @brief Returns the length of the box sides.
	 *
	 * @return
	 */
	virtual LR_vector box_sides() const = 0;

	/**
	 * @brief Returns the box's total volume.
	 *
	 * @return Simulation box's total volume
	 */
	virtual number V() = 0;

	/**
	 * @brief Returns the "absolute" position of a particle
	 *
	 * @param p pointer to the particle object
	 */
	virtual LR_vector get_abs_pos(BaseParticle *p) = 0;

	/**
	 * @brief Shifts the particle's position and stores internally the shift. Used in the fix_diffusion procedure
	 *
	 * @param p particle 
	 * @param amount displacement 
	 */
	virtual void shift_particle(BaseParticle *p, LR_vector &amount) = 0;

	const std::string INIT_EVENT = "box_initialised";
	const std::string UPDATE_EVENT = "box_updated";
};

using BoxPtr = std::shared_ptr<BaseBox>;

#endif /* BASEBOX_H_ */

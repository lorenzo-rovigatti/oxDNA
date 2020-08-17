/*
 * BaseBox.cpp
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#include "BaseBox.h"

#include "../Particles/BaseParticle.h"

BaseBox::BaseBox() {

}

BaseBox::~BaseBox() {

}

LR_vector BaseBox::min_image(const BaseParticle &p, const BaseParticle &q) {
	return min_image(p.pos, q.pos);
}

LR_vector BaseBox::min_image(const BaseParticle *p, const BaseParticle *q) {
	return min_image(p->pos, q->pos);
}

number BaseBox::sqr_min_image_distance(const BaseParticle &p, const BaseParticle &q) {
	return sqr_min_image_distance(p.pos, q.pos);
}

number BaseBox::sqr_min_image_distance(const BaseParticle *p, const BaseParticle *q) {
	return sqr_min_image_distance(p->pos, q->pos);
}

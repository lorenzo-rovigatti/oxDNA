/*
 * DNANucleotide.h
 *
 *  Created on: 04/mar/2013
 *      Author: lorenzo
 */

#ifndef DNANUCLEOTIDE_H_
#define DNANUCLEOTIDE_H_

#include "BaseParticle.h"

/**
 * @brief Represents a DNA nucleotide. Used by DNAInteraction.
 */

class DNANucleotide: public BaseParticle {
protected:
	int _grooving ;
public:
	const static LR_vector principal_axis;
	const static LR_vector stack_axis;
	const static LR_vector third_axis;

	float _POS_BACK ;
	float _POS_MM_BACK1 ;
	float _POS_MM_BACK2 ;
	float _POS_STACK ;
	float _POS_BASE ;

	enum Site {
		BACK = 0,
		STACK = 1,
		BASE = 2
	};

	DNANucleotide() = delete;
	DNANucleotide(int grooving);
	virtual ~DNANucleotide();

	virtual bool is_bonded(BaseParticle *q);
	virtual void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}
};

#endif /* DNANUCLEOTIDE_H_ */

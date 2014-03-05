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
template<typename number>
class DNANucleotide: public BaseParticle<number> {
protected:
	bool _grooving;
	LR_vector<number> _principal_axis;
	LR_vector<number> _stack_axis;
	LR_vector<number> _third_axis;
public:
	enum {
		BACK = 0,
		STACK = 1,
		BASE = 2
	};

	DNANucleotide(bool grooving);
	virtual ~DNANucleotide();

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}
};

#endif /* DNANUCLEOTIDE_H_ */

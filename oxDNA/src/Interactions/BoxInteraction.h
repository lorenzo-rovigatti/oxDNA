/*
 * BoxInteraction.h
 *
 *  Created on: 28/Oct/2013
 *      Author: Flavio
 */

#ifndef BOXINTERACTION_H_
#define BOXINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Interaction class to simulate hard boxes
 *
 *
 * The algorithm for overlap detection uses the separating axis theorem. 
 * the version present here can be probably further optimized
 * This interaction is used with interaction_type = Box
 *
 * Input options:
 *
@verbatim
box_sides = <float>, <float>, <float> (sides of the box)
@endverbatim
 */
template <typename number>
class BoxInteraction: public BaseInteraction<number, BoxInteraction<number> > {
protected:
	inline number _box_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	number _sides[3];
	number _lx, _ly, _lz;
	number _smallest, _largest;

public:
	enum {
		Box = 0
	};

	BoxInteraction();
	virtual ~BoxInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	bool generate_random_configuration_overlap (BaseParticle<number> * p, BaseParticle<number> *q, number box_side);
};

template<typename number>
number BoxInteraction<number>::_box_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	number rnorm = r->norm();
	
	// outside of the bounding sphere
	if (rnorm > this->_sqr_rcut) return (number) 0.;
	
	// inside of the inner bounding sphere
	if (rnorm < _smallest * _smallest) {
		this->set_is_infinite (true);
		return (number) 1.e12;
	}
	
	// more difficult case of checking all the 15 different cases; 
	// algorithm sketched in supp materials of 
	// Frank Smallenburg and Laura Filion's PNAS on cubes
	// http://www.pnas.org/content/suppl/2012/09/05/1211784109.DCSupplemental/pnas.1211784109_SI.pdf#STXT
	
	// here we compute the 15 potential separating axes
	LR_vector<number> sep[15];
	
	sep[0] = p->orientation.v1;
	sep[1] = p->orientation.v2;
	sep[2] = p->orientation.v3;
	
	sep[3] = q->orientation.v1;
	sep[4] = q->orientation.v2;
	sep[5] = q->orientation.v3;

	sep[6] =  p->orientation.v1.cross(q->orientation.v1);
	sep[7] =  p->orientation.v1.cross(q->orientation.v2);
	sep[8] =  p->orientation.v1.cross(q->orientation.v3);
	sep[9] =  p->orientation.v2.cross(q->orientation.v1);
	sep[10] = p->orientation.v2.cross(q->orientation.v2);
	sep[11] = p->orientation.v2.cross(q->orientation.v3);
	sep[12] = p->orientation.v3.cross(q->orientation.v1);
	sep[13] = p->orientation.v3.cross(q->orientation.v2);
	sep[14] = p->orientation.v3.cross(q->orientation.v3);

	for (int k = 6; k < 15; k ++) sep[k].normalize();

	// now we have the separating vectors; we should look for Ra and Rb
	number Ra, Rb;
	for (int k = 0; k < 15; k ++) {
		Ra = fabs(p->orientation.v1 * sep[k]) * _lx / 2.+ 
		     fabs(p->orientation.v2 * sep[k]) * _ly / 2.+
		     fabs(p->orientation.v3 * sep[k]) * _lz / 2.;

		Rb = fabs(q->orientation.v1 * sep[k]) * _lx / 2.+ 
		     fabs(q->orientation.v2 * sep[k]) * _ly / 2.+
		     fabs(q->orientation.v3 * sep[k]) * _lz / 2.;
		if (fabs((*r) * sep[k]) > (Ra + Rb)) {
			// no overlap
			return (number) 0.;
		}
	}

	// if we ended up here, it means we could not find a separating axis.	
	// this means the two particles overlap.
	this->set_is_infinite (true);
	return (number) 1.e12;
}


#endif /* BOXINTERACTION_H_ */


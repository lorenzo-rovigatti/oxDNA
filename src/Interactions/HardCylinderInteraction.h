/*
 * HardCylinderInteraction.h
 *
 *  Created on: 31/Oct/2013
 *      Author: Flavio
 */

#ifndef HARDCYLINDER_H_
#define HARDCYLINDER_H_

#include "BaseInteraction.h"
#include "InteractionUtils.h"

/**
 * @brief Interaction class to simulate Hard Cylinders.
 *
 * The cylinders are treated the intersection of a box and a spherocyliner. The
 * radius of the cylinders is assumed to be 0.5 (diameter 1.) and the length
 * (height) is set in the input file. 
 * 
 * interaction_type = HardCylinder
 *
 * @verbatim
height = <float> (cylinder length)
@endverbatim
 */

class HardCylinderInteraction: public BaseInteraction {
protected:
	inline number _hc_pot (BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _height;

public:
	enum {
		HardCylinder = 0
	};

	HardCylinderInteraction();
	virtual ~HardCylinderInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};


number HardCylinderInteraction::_hc_pot (BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	number rnorm = _computed_r.norm();
	if(rnorm > _sqr_rcut) {
		return (number) 0.f;
	}
	
	bool res1 = InteractionUtils::cylinder_overlap(p, q, _computed_r, _height);
	
	if(res1 == false) return 0.;
	set_is_infinite(true);
	return 1.e12;
}

#endif /* HARDCYLINDER_H_ */


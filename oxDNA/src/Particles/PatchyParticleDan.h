/*
 * PatchyParticleDan.h
 *
 *  Created on: 01/feb/2016
 *      Author: lorenzo -> dan
 */

#ifndef PATCHYPARTICLEDAN_H_
#define PATCHYPARTICLEDAN_H_

#include "BaseParticle.h"

/**
 * @brief Encapsulates a patchy particle. Used by PatchyInteractionDan.
 */

template<typename number>
class PatchyParticleDan : public BaseParticle<number> {
protected:
        //Base (i.e. for particles in their original orientation) patch vectors and reference vectors
	LR_vector<number> * _base_patch_vectors;
	LR_vector<number> * _base_ref_vectors;

	//bool tor_flag;

public:
        LR_vector<number> * _ref_vectors;   

	PatchyParticleDan(int N_patches, LR_vector<number> * _inp_patch_vectors, LR_vector<number> * _inp_ref_vectors, bool tor_flag);

	virtual ~PatchyParticleDan();

	void set_positions();

	virtual void copy_from(const PatchyParticleDan<number> &);

	virtual bool is_rigid_body() {
		return true;
	}

};

#endif /* PATCHYPARTICLEDAN_H_ */

/*
 * RNANucleotide.h
 *
 *  Created on: 04/mar/2013
 *      Author: petr
 */

#ifndef RNANUCLEOTIDE_H_
#define RNANUCLEOTIDE_H_

#include "BaseParticle.h"

#include "../Interactions/rna_model.h"
/**
 * @brief Represents a RNA nucleotide. Used by RNAInteraction.
 */
template<typename number>
class RNANucleotide: public BaseParticle<number> {
protected:
	static Model *model; //structure with all model constants
	LR_vector<number> _principal_axis;
	LR_vector<number> _stack_axis;
	LR_vector<number> _third_axis;

	void _set_back_position() {
		this->orientationT = this->orientation.get_transpose();

		LR_vector<number> a1 = this->orientationT.v1; // * _principal_axis;
		LR_vector<number> a3 = this->orientationT.v3; // * _stack_axis;
		LR_vector<number> a2 = this->orientationT.v2; //third_axis

		this->int_centers[BACK] =  a1 * model->RNA_POS_BACK_a1 + a2 * model->RNA_POS_BACK_a2 + a3 * model->RNA_POS_BACK_a3;
		this->int_centers[BBVECTOR_3]  = a1 * model->p3_x + a2*model->p3_y  + a3 *model->p3_z;
	    this->int_centers[BBVECTOR_5]  = a1 * model->p5_x + a2 * model->p5_y + a3 * model->p5_z;

	}
	void _set_stack_position() {
		LR_vector<number> a1 = this->orientationT.v1; // * _principal_axis;
		LR_vector<number> a2 = this->orientationT.v2; // * _stack_axis;

		/*
		a1 = orientation *_principal_axis;
		a2 = orientation *_third_axis;
		a3 = orientation *_stack_axis;
		 */

		this->int_centers[STACK_3] =  a1 * model->RNA_POS_STACK_3_a1  + a2 * model->RNA_POS_STACK_3_a2 ;
		this->int_centers[STACK_5] =  a1 * model->RNA_POS_STACK_5_a1  + a2 * model->RNA_POS_STACK_5_a2 ;

	    this->int_centers[STACK] = this->orientation*_principal_axis * (model->RNA_POS_STACK);
	}

	void _set_base_position() {this->int_centers[BASE] = this->orientation*_principal_axis* (model->RNA_POS_BASE);}
public:
	enum {
		BACK = 0,
		STACK = 1,
		BASE = 2,
		STACK_3 = 3,
		STACK_5 = 4,
		BBVECTOR_3 = 5,
		BBVECTOR_5 = 6
	};

	RNANucleotide();

	static void set_model(Model *rnamodel)
	{
			model = rnamodel;
	}

	virtual ~RNANucleotide();

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}
};

#endif /* RNANUCLEOTIDE_H_ */

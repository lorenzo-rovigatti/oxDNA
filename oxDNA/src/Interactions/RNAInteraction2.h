/**
 * @brief Handles interactions between RNA nucleotides. Contains additionally interactions (mismatch repulsion, salt) with respect to RNA Interaction class
 * Implements RNA2 model
 *
 * petr (April 2014)
 * If you want to use the RNA2 model, you need to set
 *  interaction_type = RNA2 in input file
 * Input options:
 *
 * @verbatim
 [use_average_seq = <boolean> (defaults to yes)]
 [seq_dep_file = <string> (sets the location of the files with sequence-dependent parameters)]
 [external_model = <string> (overrides default constants for the model, set in rna_model.h), by values specified by this option)]
 [salt = <float>  (sets the salt concentration in M, defaults to 1)]
 [mismatch_repulsion = <boolean> (defaults to no)]
 [mismatch_repulsion_strength = <float> (defaults to 1, sets the strength of repulsion if mismatch_repulsion is true)]
 @endverbatim
 */

#ifndef RNA2_INTERACTION_H
#define RNA2_INTERACTION_H

#include "RNAInteraction.h"
//#include "rna2_model.h"

template<typename number>
class RNA2Interaction: public RNAInteraction<number> {

private:
	number test_huckel(number rbackmod);
protected:
	//parameters of the interaction
	float _salt_concentration;
	bool _mismatch_repulsion;
	bool _debye_huckel_half_charged_ends;
	number _debye_huckel_prefactor; // this is the strength of the interaction
	number _debye_huckel_lambdafactor; //Lambda is _debye_huckel_LAMBDAFACTOR / salt^0.5
	number _debye_huckel_cutoff_factor;

	//the following values are calculated
	number _debye_huckel_Vrc; // is equal to DH(2*lambda,lambda,prefactor)
    number _debye_huckel_RC; //this is the maximum interaction distance between backbones to interact with DH
	number _debye_huckel_B; //prefactor of the quadratic cut-off
	number _debye_huckel_RHIGH; //distance after which the potential is replaced by a quadratic cut-off
	number _minus_kappa; //= -1/lambda

	//this is for the mismatch repulsion potential
	float _RNA_HYDR_MIS;
	number _fX(number r, int type,int n3, int n5);
	number _fXD(number r, int type,int n3, int n5);
	number _hydrogen_bonding_repulsion(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);


    virtual number _debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);


public:
	enum {
	  DEBYE_HUCKEL = 7
	};
	RNA2Interaction();    // Constructor
	virtual ~RNA2Interaction() {} // Destructor

	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number _hydrogen_bonding(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	virtual void get_settings(input_file &inp); //get settings from input file
	virtual void init(); // initialisation


};

#endif

#ifndef DNA3_INTERACTION_H
#define DNA3_INTERACTION_H

#include "DNA2Interaction.h"


class DNA3Interaction: public DNA2Interaction {

protected:

	/// here we store the r0, which is different in oxDNA and oxDNA2, and is set in the constructor
	number _fene_r0_SD[5][5][5][5];
	number _fene_delta_SD[5][5][5][5];
	number _fene_delta2_SD[5][5][5][5];
	number _mbf_xmax_SD[5][5][5][5];

	Mesh _mesh_f4_SD[21][5][5][5][5];
	virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	
	number _f1_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f1D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f2_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f2D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f4_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f4D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f4Dsin_SD(number t, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f5_SD(number f, int type, int n3, int n5);
	number _f5D_SD(number f, int type, int n3, int n5);
	
	number _fakef4_SD(number t, void * par);
	number _fakef4D_SD(number t, void * par);
	
	//virtual number _custom_f4_SD (number cost, int i, int n3, int n5) { return this->_mesh_f4_SD[i][n3][n5].query(cost); }
	virtual number _custom_f4_SD (number cost, int i, int n3_2, int n3_1, int n5_1, int n5_2) { return this->_mesh_f4_SD[i][n3_2][n3_1][n5_1][n5_2].query(cost); }

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	//virtual number _custom_f4D_SD (number cost, int i, int n3, int n5) { return this->_mesh_f4_SD[i][n3][n5].query_derivative(cost); }
	virtual number _custom_f4D_SD (number cost, int i, int n3_2, int n3_1, int n5_1, int n5_2) { return this->_mesh_f4_SD[i][n3_2][n3_1][n5_1][n5_2].query_derivative(cost); }

public:
	DNA3Interaction();
	
	virtual void init();
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
	
	number F1_SD_A[2][5][5][5][5];
	number F1_SD_RC[2][5][5][5][5];
	number F1_SD_R0[2][5][5][5][5];
	number F1_SD_BLOW[2][5][5][5][5];
	number F1_SD_BHIGH[2][5][5][5][5];
	number F1_SD_RLOW[2][5][5][5][5];
	number F1_SD_RHIGH[2][5][5][5][5];
	number F1_SD_RCLOW[2][5][5][5][5];
	number F1_SD_RCHIGH[2][5][5][5][5];
	number F1_SD_SHIFT[2][5][5][5][5];

	number F2_SD_K[4][5][5][5][5];
	number F2_SD_RC[4][5][5][5][5];
	number F2_SD_R0[4][5][5][5][5];
	number F2_SD_BLOW[4][5][5][5][5];
	number F2_SD_RLOW[4][5][5][5][5];
	number F2_SD_RCLOW[4][5][5][5][5];
	number F2_SD_BHIGH[4][5][5][5][5];
	number F2_SD_RCHIGH[4][5][5][5][5];
	number F2_SD_RHIGH[4][5][5][5][5];

	number F4_SD_THETA_A[21][5][5][5][5];
	number F4_SD_THETA_B[21][5][5][5][5];
	number F4_SD_THETA_T0[21][5][5][5][5];
	number F4_SD_THETA_TS[21][5][5][5][5];
	number F4_SD_THETA_TC[21][5][5][5][5];

	number F5_SD_PHI_A[4][5][5];
	number F5_SD_PHI_B[4][5][5];
	number F5_SD_PHI_XC[4][5][5];
	number F5_SD_PHI_XS[4][5][5];
	
	number MESH_F4_SD_POINTS[21][5][5][5][5];	
	
	
        ~DNA3Interaction();
};

#endif /* DNA3_INTERACTION_H */

#ifndef DNA3_INTERACTION_H
#define DNA3_INTERACTION_H

#define TETRAMER_DIM_A 6 // this is the number of nucleotide types including dummies
#define TETRAMER_DIM_B 5 // this is the number of nucleotide types

#define NO_TYPE 5

#include "DNA2Interaction.h"
#include "../Utilities/oxdna3_utils.h"

class DNA3Interaction: public DNA2Interaction {

protected:

	/// here we store the r0, which is different in oxDNA and oxDNA2, and is set in the constructor
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _fene_r0_SD;
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _fene_delta_SD;
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _fene_delta2_SD;
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _mbf_xmax_SD;
    number _fene_eps;

	Mesh _mesh_f4_SD[21][TETRAMER_DIM_A][TETRAMER_DIM_B][TETRAMER_DIM_B][TETRAMER_DIM_A]; // this is not used by the GPU so we can keep using C-style arrays
	virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _f1_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f1D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f2_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f2D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f4_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f4D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f4Dsin_SD(number t, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f5_SD(number f, int type, int n3_2, int n3_1, int n5_1, int n5_2);
	number _f5D_SD(number f, int type, int n3_2, int n3_1, int n5_1, int n5_2);

	number _fakef4_SD(number t, void * par);
	number _fakef4D_SD(number t, void * par);

	virtual number _custom_f4_SD(number cost, int i, int n3_2, int n3_1, int n5_1, int n5_2) { return _mesh_f4_SD[i][n3_2][n3_1][n5_1][n5_2].query(cost); }

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4D_SD(number cost, int i, int n3_2, int n3_1, int n5_1, int n5_2) { return _mesh_f4_SD[i][n3_2][n3_1][n5_1][n5_2].query_derivative(cost); }

public:
	DNA3Interaction();

	virtual void init();
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_A[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_RC[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_R0[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_BLOW[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_BHIGH[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_RLOW[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_RHIGH[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_RCLOW[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_RCHIGH[2];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F1_SD_SHIFT[2];

	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_K[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_RC[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_R0[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_BLOW[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_RLOW[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_RCLOW[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_BHIGH[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_RCHIGH[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F2_SD_RHIGH[4];

	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F4_SD_THETA_A[21];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F4_SD_THETA_B[21];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F4_SD_THETA_T0[21];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F4_SD_THETA_TS[21];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F4_SD_THETA_TC[21];

	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F5_SD_PHI_A[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F5_SD_PHI_B[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F5_SD_PHI_XC[4];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> F5_SD_PHI_XS[4];

    MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _excl_s[7];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _excl_r[7];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _excl_b[7];
	MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> _excl_rc[7];

	std::array<int, 21> MESH_F4_SD_POINTS;

    number get_fene_delta2_SD(int ty1, int ty2, int ty3, int ty4) { return _fene_delta2_SD(ty1, ty2, ty3, ty4); }
    number get_fene_delta_SD(int ty1, int ty2, int ty3, int ty4) { return _fene_delta_SD(ty1, ty2, ty3, ty4); }
    number get_fene_r0_SD(int ty1, int ty2, int ty3, int ty4)  { return _fene_r0_SD(ty1, ty2, ty3, ty4); }

    ~DNA3Interaction();
};

/**
 * @brief Handles interactions between DNA3 nucleotides without using meshes.
 *
 * This interaction is selected with
 * interaction_type = DNA3_nomesh
 */
class DNA3Interaction_nomesh: public DNA3Interaction {
protected:

	number _custom_f4_SD(number cost, int i, int n3_2, int n3_1, int n5_1, int n5_2) override {
		return this->_fakef4_SD(cost, std::array<int, 5>({i, n3_2, n3_1, n5_1, n5_2}).data());
	}

	number _custom_f4D_SD(number cost, int i, int n3_2, int n3_1, int n5_1, int n5_2) override {
		return this->_fakef4D_SD(cost, std::array<int, 5>({i, n3_2, n3_1, n5_1, n5_2}).data());
	}

public:
	DNA3Interaction_nomesh() {
	}
	
	virtual ~DNA3Interaction_nomesh() {
	}
};

#endif /* DNA3_INTERACTION_H */

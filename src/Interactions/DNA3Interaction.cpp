#include "DNA3Interaction.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../Particles/DNANucleotide.h"

DNA3Interaction::DNA3Interaction() : DNA2Interaction() {
    OX_LOG(Logger::LOG_INFO, "Using the experimental DNA3Interaction");
    _grooving = 2;

    _fene_eps = FENE_EPS;

    MESH_F4_SD_POINTS[HYDR_F4_THETA1] = HYDR_T1_MESH_POINTS;
    MESH_F4_SD_POINTS[HYDR_F4_THETA2] = HYDR_T2_MESH_POINTS;
    MESH_F4_SD_POINTS[HYDR_F4_THETA4] = HYDR_T4_MESH_POINTS;
    MESH_F4_SD_POINTS[HYDR_F4_THETA7] = HYDR_T7_MESH_POINTS;

    MESH_F4_SD_POINTS[STCK_F4_THETA4] = STCK_T4_MESH_POINTS;
    MESH_F4_SD_POINTS[STCK_F4_THETA5] = STCK_T5_MESH_POINTS;

    MESH_F4_SD_POINTS[CRST_F4_THETA1] = CRST_T1_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA2] = CRST_T2_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA4] = CRST_T4_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA7] = CRST_T7_MESH_POINTS;

    MESH_F4_SD_POINTS[CXST_F4_THETA1] = CXST_T1_MESH_POINTS;
    MESH_F4_SD_POINTS[CXST_F4_THETA4] = CXST_T4_MESH_POINTS;
    MESH_F4_SD_POINTS[CXST_F4_THETA5] = CXST_T5_MESH_POINTS;

    MESH_F4_SD_POINTS[CRST_F4_THETA1_33] = CRST_T1_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA2_33] = CRST_T2_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA4_33] = CRST_T4_33_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA7_33] = CRST_T7_MESH_POINTS;

    MESH_F4_SD_POINTS[CRST_F4_THETA1_55] = CRST_T1_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA2_55] = CRST_T2_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA4_55] = CRST_T4_55_MESH_POINTS;
    MESH_F4_SD_POINTS[CRST_F4_THETA7_55] = CRST_T7_MESH_POINTS;

    // initialise all sequence dependent parameters with their average value.
    // these values are overwritten with the sequence dependent values in init().
    _fene_r0_SD.fill(FENE_R0_OXDNA2);
    _fene_delta_SD.fill(FENE_DELTA);
    _fene_delta2_SD.fill(FENE_DELTA2);

    _excl_s[0].fill(EXCL_S1);
    _excl_s[1].fill(EXCL_S2);
    _excl_s[2].fill(EXCL_S3);
    _excl_s[3].fill(EXCL_S4);
    _excl_s[4].fill(EXCL_S5);
    _excl_s[5].fill(EXCL_S6);
    _excl_s[6].fill(EXCL_S7);

    _excl_r[0].fill(EXCL_R1);
    _excl_r[1].fill(EXCL_R2);
    _excl_r[2].fill(EXCL_R3);
    _excl_r[3].fill(EXCL_R4);
    _excl_r[4].fill(EXCL_R5);
    _excl_r[5].fill(EXCL_R6);
    _excl_r[6].fill(EXCL_R7);

    _excl_b[0].fill(EXCL_B1);
    _excl_b[1].fill(EXCL_B2);
    _excl_b[2].fill(EXCL_B3);
    _excl_b[3].fill(EXCL_B4);
    _excl_b[4].fill(EXCL_B5);
    _excl_b[5].fill(EXCL_B6);
    _excl_b[6].fill(EXCL_B7);

    _excl_rc[0].fill(EXCL_RC1);
    _excl_rc[1].fill(EXCL_RC2);
    _excl_rc[2].fill(EXCL_RC3);
    _excl_rc[3].fill(EXCL_RC4);
    _excl_rc[4].fill(EXCL_RC5);
    _excl_rc[5].fill(EXCL_RC6);
    _excl_rc[6].fill(EXCL_RC7);


    F1_SD_EPS[0].fill(HYDR_EPS_OXDNA);
    F1_SD_EPS[1].fill(STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * _T);

    F1_SD_A[0].fill(HYDR_A);
    F1_SD_A[1].fill(STCK_A);

    F1_SD_RC[0].fill(HYDR_RC);
    F1_SD_RC[1].fill(STCK_RC);

    F1_SD_R0[0].fill(HYDR_R0);
    F1_SD_R0[1].fill(STCK_R0);

    F1_SD_BLOW[0].fill(HYDR_BLOW);
    F1_SD_BLOW[1].fill(STCK_BLOW);

    F1_SD_BHIGH[0].fill(HYDR_BHIGH);
    F1_SD_BHIGH[1].fill(STCK_BHIGH);

    F1_SD_RLOW[0].fill(HYDR_RLOW);
    F1_SD_RLOW[1].fill(STCK_RLOW);

    F1_SD_RHIGH[0].fill(HYDR_RHIGH);
    F1_SD_RHIGH[1].fill(STCK_RHIGH);

    F1_SD_RCLOW[0].fill(HYDR_RCLOW);
    F1_SD_RCLOW[1].fill(STCK_RCLOW);

    F1_SD_RCHIGH[0].fill(HYDR_RCHIGH);
    F1_SD_RCHIGH[1].fill(STCK_RCHIGH);

    F4_SD_THETA_A[0].fill(STCK_THETA4_A);
    F4_SD_THETA_A[1].fill(STCK_THETA5_A);

    F4_SD_THETA_A[0].fill(STCK_THETA4_A);
    F4_SD_THETA_A[1].fill(STCK_THETA5_A);

    F4_SD_THETA_A[2].fill(HYDR_THETA1_A);
    F4_SD_THETA_A[3].fill(HYDR_THETA2_A);
    F4_SD_THETA_A[4].fill(HYDR_THETA4_A);
    F4_SD_THETA_A[5].fill(HYDR_THETA7_A);

    F4_SD_THETA_A[6].fill(CRST_THETA1_A);
    F4_SD_THETA_A[7].fill(CRST_THETA2_A);
    F4_SD_THETA_A[8].fill(CRST_THETA4_A);
    F4_SD_THETA_A[9].fill(CRST_THETA7_A);

    F4_SD_THETA_A[10].fill(CXST_THETA1_A);
    F4_SD_THETA_A[11].fill(CXST_THETA4_A);
    F4_SD_THETA_A[12].fill(CXST_THETA5_A);

    F4_SD_THETA_A[13].fill(CRST_THETA1_A_33);
    F4_SD_THETA_A[14].fill(CRST_THETA2_A_33);
    F4_SD_THETA_A[15].fill(CRST_THETA4_A_33);
    F4_SD_THETA_A[16].fill(CRST_THETA7_A_33);
    F4_SD_THETA_A[17].fill(CRST_THETA1_A_55);
    F4_SD_THETA_A[18].fill(CRST_THETA2_A_55);
    F4_SD_THETA_A[19].fill(CRST_THETA4_A_55);
    F4_SD_THETA_A[20].fill(CRST_THETA7_A_55);

    F4_SD_THETA_B[0].fill(STCK_THETA4_B);
    F4_SD_THETA_B[1].fill(STCK_THETA5_B);

    F4_SD_THETA_B[2].fill(HYDR_THETA1_B);
    F4_SD_THETA_B[3].fill(HYDR_THETA2_B);
    F4_SD_THETA_B[4].fill(HYDR_THETA4_B);
    F4_SD_THETA_B[5].fill(HYDR_THETA7_B);

    F4_SD_THETA_B[6].fill(CRST_THETA1_B);
    F4_SD_THETA_B[7].fill(CRST_THETA2_B);
    F4_SD_THETA_B[8].fill(CRST_THETA4_B);
    F4_SD_THETA_B[9].fill(CRST_THETA7_B);

    F4_SD_THETA_B[10].fill(CXST_THETA1_B);
    F4_SD_THETA_B[11].fill(CXST_THETA4_B);
    F4_SD_THETA_B[12].fill(CXST_THETA5_B);

    F4_SD_THETA_B[13].fill(CRST_THETA1_B_33);
    F4_SD_THETA_B[14].fill(CRST_THETA2_B_33);
    F4_SD_THETA_B[15].fill(CRST_THETA4_B_33);
    F4_SD_THETA_B[16].fill(CRST_THETA7_B_33);
    F4_SD_THETA_B[17].fill(CRST_THETA1_B_55);
    F4_SD_THETA_B[18].fill(CRST_THETA2_B_55);
    F4_SD_THETA_B[19].fill(CRST_THETA4_B_55);
    F4_SD_THETA_B[20].fill(CRST_THETA7_B_55);

    F4_SD_THETA_T0[0].fill(STCK_THETA4_T0);
    F4_SD_THETA_T0[1].fill(STCK_THETA5_T0);

    F4_SD_THETA_T0[2].fill(HYDR_THETA1_T0);
    F4_SD_THETA_T0[3].fill(HYDR_THETA2_T0);
    F4_SD_THETA_T0[4].fill(HYDR_THETA4_T0);
    F4_SD_THETA_T0[5].fill(HYDR_THETA7_T0);

    F4_SD_THETA_T0[6].fill(CRST_THETA1_T0);
    F4_SD_THETA_T0[7].fill(CRST_THETA2_T0);
    F4_SD_THETA_T0[8].fill(CRST_THETA4_T0);
    F4_SD_THETA_T0[9].fill(CRST_THETA7_T0);

    F4_SD_THETA_T0[10].fill(CXST_THETA1_T0_OXDNA);
    F4_SD_THETA_T0[11].fill(CXST_THETA4_T0);
    F4_SD_THETA_T0[12].fill(CXST_THETA5_T0);
    F4_SD_THETA_T0[13].fill(CRST_THETA1_T0_33);
    F4_SD_THETA_T0[14].fill(CRST_THETA2_T0_33);
    F4_SD_THETA_T0[15].fill(CRST_THETA4_T0_33);
    F4_SD_THETA_T0[16].fill(CRST_THETA7_T0_33);
    F4_SD_THETA_T0[17].fill(CRST_THETA1_T0_55);
    F4_SD_THETA_T0[18].fill(CRST_THETA2_T0_55);
    F4_SD_THETA_T0[19].fill(CRST_THETA4_T0_55);
    F4_SD_THETA_T0[20].fill(CRST_THETA7_T0_55);

    F4_SD_THETA_TS[0].fill(STCK_THETA4_TS);
    F4_SD_THETA_TS[1].fill(STCK_THETA5_TS);

    F4_SD_THETA_TS[2].fill(HYDR_THETA1_TS);
    F4_SD_THETA_TS[3].fill(HYDR_THETA2_TS);
    F4_SD_THETA_TS[4].fill(HYDR_THETA4_TS);
    F4_SD_THETA_TS[5].fill(HYDR_THETA7_TS);

    F4_SD_THETA_TS[6].fill(CRST_THETA1_TS);
    F4_SD_THETA_TS[7].fill(CRST_THETA2_TS);
    F4_SD_THETA_TS[8].fill(CRST_THETA4_TS);
    F4_SD_THETA_TS[9].fill(CRST_THETA7_TS);

    F4_SD_THETA_TS[10].fill(CXST_THETA1_TS);
    F4_SD_THETA_TS[11].fill(CXST_THETA4_TS);
    F4_SD_THETA_TS[12].fill(CXST_THETA5_TS);

    F4_SD_THETA_TS[13].fill(CRST_THETA1_TS_33);
    F4_SD_THETA_TS[14].fill(CRST_THETA2_TS_33);
    F4_SD_THETA_TS[15].fill(CRST_THETA4_TS_33);
    F4_SD_THETA_TS[16].fill(CRST_THETA7_TS_33);
    F4_SD_THETA_TS[17].fill(CRST_THETA1_TS_55);
    F4_SD_THETA_TS[18].fill(CRST_THETA2_TS_55);
    F4_SD_THETA_TS[19].fill(CRST_THETA4_TS_55);
    F4_SD_THETA_TS[20].fill(CRST_THETA7_TS_55);

    F4_SD_THETA_TC[0].fill(STCK_THETA4_TC);
    F4_SD_THETA_TC[1].fill(STCK_THETA5_TC);

    F4_SD_THETA_TC[2].fill(HYDR_THETA1_TC);
    F4_SD_THETA_TC[3].fill(HYDR_THETA2_TC);
    F4_SD_THETA_TC[4].fill(HYDR_THETA4_TC);
    F4_SD_THETA_TC[5].fill(HYDR_THETA7_TC);

    F4_SD_THETA_TC[6].fill(CRST_THETA1_TC);
    F4_SD_THETA_TC[7].fill(CRST_THETA2_TC);
    F4_SD_THETA_TC[8].fill(CRST_THETA4_TC);
    F4_SD_THETA_TC[9].fill(CRST_THETA7_TC);

    F4_SD_THETA_TC[10].fill(CXST_THETA1_TC);
    F4_SD_THETA_TC[11].fill(CXST_THETA4_TC);
    F4_SD_THETA_TC[12].fill(CXST_THETA5_TC);

    F4_SD_THETA_TC[13].fill(CRST_THETA1_TC_33);
    F4_SD_THETA_TC[14].fill(CRST_THETA2_TC_33);
    F4_SD_THETA_TC[15].fill(CRST_THETA4_TC_33);
    F4_SD_THETA_TC[16].fill(CRST_THETA7_TC_33);
    F4_SD_THETA_TC[17].fill(CRST_THETA1_TC_55);
    F4_SD_THETA_TC[18].fill(CRST_THETA2_TC_55);
    F4_SD_THETA_TC[19].fill(CRST_THETA4_TC_55);
    F4_SD_THETA_TC[20].fill(CRST_THETA7_TC_55);

    F2_SD_K[0].fill(CRST_K);
    F2_SD_K[1].fill(CXST_K_OXDNA);
    F2_SD_K[2].fill(CRST_K_33);
    F2_SD_K[3].fill(CRST_K_55);

    F2_SD_RC[0].fill(CRST_RC);
    F2_SD_RC[1].fill(CXST_RC);
    F2_SD_RC[2].fill(CRST_RC_33);
    F2_SD_RC[3].fill(CRST_RC_55);

    F2_SD_R0[0].fill(CRST_R0);
    F2_SD_R0[1].fill(CXST_R0);
    F2_SD_R0[2].fill(CRST_R0_33);
    F2_SD_R0[3].fill(CRST_R0_55);

    F2_SD_BLOW[0].fill(CRST_BLOW);
    F2_SD_BLOW[1].fill(CXST_BLOW);
    F2_SD_BLOW[2].fill(CRST_BLOW_33);
    F2_SD_BLOW[3].fill(CRST_BLOW_55);

    F2_SD_BHIGH[0].fill(CRST_BHIGH);
    F2_SD_BHIGH[1].fill(CXST_BHIGH);
    F2_SD_BHIGH[2].fill(CRST_BHIGH_33);
    F2_SD_BHIGH[3].fill(CRST_BHIGH_55);

    F2_SD_RLOW[0].fill(CRST_RLOW);
    F2_SD_RLOW[1].fill(CXST_RLOW);
    F2_SD_RLOW[2].fill(CRST_RLOW_33);
    F2_SD_RLOW[3].fill(CRST_RLOW_55);

    F2_SD_RHIGH[0].fill(CRST_RHIGH);
    F2_SD_RHIGH[1].fill(CXST_RHIGH);
    F2_SD_RHIGH[2].fill(CRST_RHIGH_33);
    F2_SD_RHIGH[3].fill(CRST_RHIGH_55);

    F2_SD_RCLOW[0].fill(CRST_RCLOW);
    F2_SD_RCLOW[1].fill(CXST_RCLOW);
    F2_SD_RCLOW[2].fill(CRST_RCLOW_33);
    F2_SD_RCLOW[3].fill(CRST_RCLOW_55);

    F2_SD_RCHIGH[0].fill(CRST_RCHIGH);
    F2_SD_RCHIGH[1].fill(CXST_RCHIGH);
    F2_SD_RCHIGH[2].fill(CRST_RCHIGH_33);
    F2_SD_RCHIGH[3].fill(CRST_RCHIGH_55);

    F5_SD_PHI_A[0].fill(STCK_PHI1_A);
    F5_SD_PHI_A[1].fill(STCK_PHI2_A);
    F5_SD_PHI_A[2].fill(CXST_PHI3_A);
    F5_SD_PHI_A[3].fill(CXST_PHI4_A);

    F5_SD_PHI_B[0].fill(STCK_PHI1_B);
    F5_SD_PHI_B[1].fill(STCK_PHI2_B);
    F5_SD_PHI_B[2].fill(CXST_PHI3_B);
    F5_SD_PHI_B[3].fill(CXST_PHI3_B);

    F5_SD_PHI_XC[0].fill(STCK_PHI1_XC);
    F5_SD_PHI_XC[1].fill(STCK_PHI2_XC);
    F5_SD_PHI_XC[2].fill(CXST_PHI3_XC);
    F5_SD_PHI_XC[3].fill(CXST_PHI4_XC);

    F5_SD_PHI_XS[0].fill(STCK_PHI1_XS);
    F5_SD_PHI_XS[1].fill(STCK_PHI2_XS);
    F5_SD_PHI_XS[2].fill(CXST_PHI3_XS);
    F5_SD_PHI_XS[3].fill(CXST_PHI4_XS);
}

void DNA3Interaction::init() {
    Logger::instance()->disable_log("DNA3Interaction");
    DNA2Interaction::init();
	Logger::instance()->enable_log("DNA3Interaction");

    // if we want to use the sequence dependence then we overwrite all the
    // parameters regarding interactions between true bases (i.e. the
    // interactions between dummy bases or regular bases and dummy bases
    // keeps the default value)
    // the epsilons (Stacking and HB) are overwritten in DNAInteraction,
    // here we overwrite all the other parameters

    input_file seq_file;
    seq_file.init_from_filename(_seq_filename.c_str());
    char key[256];
    float tmp_value;

    float stck_fact_eps;
    getInputFloat(&seq_file, "STCK_FACT_EPS", &stck_fact_eps, 1);

    sprintf(key, "FENE_EPS");
    if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_eps = tmp_value;

    if(!_average) {
        if(seq_file.state == ERROR) {
            throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename.c_str());
        }

        // read independent parameters from SD file
        // independent = not set by continuity and differentiability
        for(int i = 0; i < TETRAMER_DIM_A - 2; i++) {
            for(int j = 0; j < TETRAMER_DIM_B - 1; j++) {
                for(int k = 0; k < TETRAMER_DIM_B - 1; k++) {
                    for(int l = 0; l < TETRAMER_DIM_A - 2; l++) {
                        // EXCLUDED VOLUME

                        sprintf(key, "EXCL_S1_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[0](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_S2_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[1](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_S3_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[2](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_S4_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[3](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_S5_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[4](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_S6_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[5](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_S7_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_s[6](i, j, k, l) = tmp_value;

                        sprintf(key, "EXCL_R1_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[0](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_R2_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[1](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_R3_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[2](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_R4_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[3](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_R5_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[4](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_R6_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[5](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_R7_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_r[6](i, j, k, l) = tmp_value;

                        /*
                        sprintf(key, "EXCL_B1_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[0](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_B2_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[1](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_B3_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[2](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_B4_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[3](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_B5_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[4](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_B6_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[5](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_B7_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_b[6](i, j, k, l) = tmp_value;

                        sprintf(key, "EXCL_RC1_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[0](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_RC2_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[1](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_RC3_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[2](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_RC4_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[3](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_RC5_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[4](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_RC6_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[5](i, j, k, l) = tmp_value;
                        sprintf(key, "EXCL_RC7_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _excl_rc[6](i, j, k, l) = tmp_value;
                        */

                        // FENE

                        sprintf(key, "FENE_R0_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_r0_SD(i, j, k, l) = tmp_value;
                        sprintf(key, "FENE_DELTA_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_delta_SD(i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "FENE_DELTA2_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) _fene_delta2_SD(i, j, k, l) = tmp_value;
                                */

                        // F1

                        // this is pair type dependent
                        sprintf(key, "HYDR_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_EPS[HYDR_F1](k, i, j, l) = tmp_value;

                        //tetramer type dependent
                        sprintf(key, "HYDR_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_A[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_RC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RC[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_R0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_R0[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_RLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RLOW[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_RHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RHIGH[HYDR_F1](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_RCLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCLOW[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_RCHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCHIGH[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_BLOW_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BLOW[HYDR_F1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_BHIGH_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BHIGH[HYDR_F1](k, i, j, l) = tmp_value;
                        */

                        // update shift
                        // F1_SD_SHIFT[HYDR_F1](i, j, k, l) = F1_EPS[HYDR_F1][j][k] * SQR(1 - exp(-(F1_SD_RC[HYDR_F1](i, j, k, l) - F1_SD_R0[HYDR_F1](i, j, k, l)) * F1_SD_A[HYDR_F1](i, j, k, l)));

                        // STACKING

                        sprintf(key, "STCK_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_EPS[STCK_F1](i, j, k, l) = tmp_value * (1.0 - stck_fact_eps + (_T * 9.0 * stck_fact_eps));


                        sprintf(key, "STCK_A_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_A[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_RC_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RC[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_R0_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_R0[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_RLOW_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RLOW[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_RHIGH_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RHIGH[STCK_F1](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "STCK_RCLOW_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCLOW[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_RCHIGH_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_RCHIGH[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_BLOW_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BLOW[STCK_F1](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_BHIGH_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j),Utils::encode_base(k),Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F1_SD_BHIGH[STCK_F1](i, j, k, l) = tmp_value;
                        */

                        // update shift
                        // F1_SD_SHIFT[STCK_F1](i, j, k, l) = F1_EPS[STCK_F1][j][k] * SQR(1 - exp(-(F1_SD_RC[STCK_F1](i, j, k, l) - F1_SD_R0[STCK_F1](i, j, k, l)) * F1_SD_A[STCK_F1](i, j, k, l)));

                        // F4

                        // HYDROGEN
                        sprintf(key, "HYDR_THETA1_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA1_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA1_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA1](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_THETA1_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA1](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA1_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA1](k, i, j, l) = tmp_value;
                        */

                        sprintf(key, "HYDR_THETA2_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA2](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA2_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA2](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA2_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA2](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_THETA2_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA2](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA2_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA2](k, i, j, l) = tmp_value;
                        */

                        sprintf(key, "HYDR_THETA3_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA3](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA3_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA3](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA3_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA3](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_THETA3_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA3](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA3_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA3](k, i, j, l) = tmp_value;
                        */

                        sprintf(key, "HYDR_THETA4_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA4](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA4_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA4](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA4_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA4](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_THETA4_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA4](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA4_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA4](k, i, j, l) = tmp_value;
                        */

                        sprintf(key, "HYDR_THETA7_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA7](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA7_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA7](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA7_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA7](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_THETA7_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA7](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA7_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA7](k, i, j, l) = tmp_value;
                        */

                        sprintf(key, "HYDR_THETA8_A_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[HYDR_F4_THETA8](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA8_T0_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[HYDR_F4_THETA8](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA8_TS_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[HYDR_F4_THETA8](k, i, j, l) = tmp_value;
                        /*
                        sprintf(key, "HYDR_THETA8_TC_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[HYDR_F4_THETA8](k, i, j, l) = tmp_value;
                        sprintf(key, "HYDR_THETA8_B_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[HYDR_F4_THETA8](k, i, j, l) = tmp_value;
                        */

                        // STACKING

                        sprintf(key, "STCK_THETA4_A_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[STCK_F4_THETA4](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_THETA4_T0_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[STCK_F4_THETA4](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_THETA4_TS_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[STCK_F4_THETA4](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "STCK_THETA4_TC_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[STCK_F4_THETA4](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_THETA4_B_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[STCK_F4_THETA4](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "STCK_THETA5_A_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[STCK_F4_THETA5](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_THETA5_T0_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[STCK_F4_THETA5](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_THETA5_TS_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[STCK_F4_THETA5](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "STCK_THETA5_TC_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[STCK_F4_THETA5](i, j, k, l) = tmp_value;
                        sprintf(key, "STCK_THETA5_B_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[STCK_F4_THETA5](i, j, k, l) = tmp_value;
                        */

                        // CROSS STACKING

                        // tetramer dependent

                        sprintf(key, "CRST_THETA4_A_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA4_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA4_T0_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA4_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA4_TS_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA4_33](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA4_TC_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA4_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA4_B_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA4_33](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA4_A_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA4_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA4_T0_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA4_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA4_TS_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA4_55](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA4_TC_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA4_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA4_B_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA4_55](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA1_A_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA1_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA1_T0_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA1_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA1_TS_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA1_33](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA1_TC_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA1_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA1_B_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA1_33](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA2_A_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA2_T0_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA2_TS_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA2_33](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA2_TC_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA2_B_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA2_33](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA7_A_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA7_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA7_T0_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA7_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA7_TS_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA7_33](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA7_TC_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA7_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA7_B_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA7_33](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA1_A_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA1_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA1_T0_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA1_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA1_TS_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA1_55](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA1_TC_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA1_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA1_B_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA1_55](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA2_A_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA2_T0_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA2_TS_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA2_55](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA2_TC_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA2_B_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA2_55](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_THETA7_A_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_A[CRST_F4_THETA7_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA7_T0_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_T0[CRST_F4_THETA7_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA7_TS_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TS[CRST_F4_THETA7_55](i, j, k, l) = tmp_value;
                        /*
                        sprintf(key, "CRST_THETA7_TC_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_TC[CRST_F4_THETA7_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_THETA7_B_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F4_SD_THETA_B[CRST_F4_THETA7_55](i, j, k, l) = tmp_value;
                        */

                        // F2

                        // CROSS_STACKING
                        // SYMMETRIC -- TODO.
                        // ASYMMETRIC

                        // this is pair type dependent
                        sprintf(key, "CRST_K_33_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_K[CRST_F2_33](k, i, j, l) = tmp_value;
                        sprintf(key, "CRST_K_55_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_K[CRST_F2_55](k, i, j, l) = tmp_value;

                        // tetramer dependent
                        sprintf(key, "CRST_R0_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_R0[CRST_F2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RC_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RC[CRST_F2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RLOW_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RLOW[CRST_F2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RHIGH_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RHIGH[CRST_F2_33](i, j, k, l) = tmp_value;

                        /*
                        sprintf(key, "CRST_BLOW_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BLOW[CRST_F2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_BHIGH_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BHIGH[CRST_F2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RCLOW_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCLOW[CRST_F2_33](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RCHIGH_33_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCHIGH[CRST_F2_33](i, j, k, l) = tmp_value;
                        */

                        sprintf(key, "CRST_R0_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_R0[CRST_F2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RC_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RC[CRST_F2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RLOW_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RLOW[CRST_F2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RHIGH_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RHIGH[CRST_F2_55](i, j, k, l) = tmp_value;

                        /*
                        sprintf(key, "CRST_BLOW_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BLOW[CRST_F2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_BHIGH_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_BHIGH[CRST_F2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RCLOW_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCLOW[CRST_F2_55](i, j, k, l) = tmp_value;
                        sprintf(key, "CRST_RCHIGH_55_%c_%c_%c_%c", Utils::encode_base(i), Utils::encode_base(j), Utils::encode_base(k), Utils::encode_base(l));
                        if(getInputFloat(&seq_file, key, &tmp_value, 0) == KEY_FOUND) F2_SD_RCHIGH[CRST_F2_55](i, j, k, l) = tmp_value;
                        */

                        // F5 -- TODO, eventually

                        // MESH POINTS --TODO, ?
                    }
                }
            }
        }

        // Set parameters for ends junctions
        for(int i = 0; i < TETRAMER_DIM_A; i++) {
            for(int j = 0; j < TETRAMER_DIM_B - 1; j++) {
                for(int k = 0; k < TETRAMER_DIM_B - 1; k++) {
                    // For 2d parameters we just copy the value for 0ij0 (par_kijl is the same for every k,l)

                    F2_SD_K[CRST_F2_33](i, j, k, TETRAMER_DIM_A - 1) = F2_SD_K[CRST_F2_33](0, j, k, 0);
                    F2_SD_K[CRST_F2_33](TETRAMER_DIM_A - 1, j, k, i) = F2_SD_K[CRST_F2_33](0, j, k, 0);

                    F2_SD_K[CRST_F2_55](i, j, k, TETRAMER_DIM_A - 1) = F2_SD_K[CRST_F2_55](0, j, k, 0);
                    F2_SD_K[CRST_F2_55](TETRAMER_DIM_A - 1, j, k, i) = F2_SD_K[CRST_F2_55](0, j, k, 0);

                    // fene
                    _fene_delta_SD(i, j, k, TETRAMER_DIM_A - 1) = _fene_delta_SD.get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                    // if(i==0 && j ==0 && k == 0) std::cout << "Fene delta ends: " << _fene_delta_SD[i][j][k][5] << std::endl;
                    _fene_delta_SD(TETRAMER_DIM_A - 1, j, k, i) = _fene_delta_SD.get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                    _fene_r0_SD(i, j, k, TETRAMER_DIM_A - 1) = _fene_r0_SD.get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                    // if(i==0 && j ==0 && k == 0) std::cout << "Fene r0 ends: " << _fene_r0_SD[i][j][k][5] << std::endl;
                    _fene_r0_SD(TETRAMER_DIM_A - 1, j, k, i) = _fene_r0_SD.get_average_par(TETRAMER_DIM_A - 1, j, k ,i);

                    // f1
                    for(int m = 0; m < 2; m++) {
                        F1_SD_EPS[m](i, j, k, TETRAMER_DIM_A - 1) = F1_SD_EPS[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F1_SD_A[m](i, j, k, TETRAMER_DIM_A - 1) = F1_SD_A[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F1_SD_R0[m](i, j, k, TETRAMER_DIM_A - 1) = F1_SD_R0[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F1_SD_RC[m](i, j, k, TETRAMER_DIM_A - 1) = F1_SD_RC[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F1_SD_RLOW[m](i, j, k, TETRAMER_DIM_A - 1) = F1_SD_RLOW[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F1_SD_RHIGH[m](i, j, k, TETRAMER_DIM_A - 1) = F1_SD_RHIGH[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);

                        F1_SD_EPS[m](TETRAMER_DIM_A - 1, j, k, i) = F1_SD_EPS[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F1_SD_A[m](TETRAMER_DIM_A - 1, j, k, i) = F1_SD_A[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F1_SD_R0[m](TETRAMER_DIM_A - 1, j, k, i) = F1_SD_R0[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F1_SD_RC[m](TETRAMER_DIM_A - 1, j, k, i) = F1_SD_RC[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F1_SD_RLOW[m](TETRAMER_DIM_A - 1, j, k, i) = F1_SD_RLOW[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F1_SD_RHIGH[m](TETRAMER_DIM_A - 1, j, k, i) = F1_SD_RHIGH[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                    }

                    // f2
                    for(int m = 0; m < 4; m++) {
                        F2_SD_R0[m](i, j, k, TETRAMER_DIM_A - 1) = F2_SD_R0[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F2_SD_RC[m](i, j, k, TETRAMER_DIM_A - 1) = F2_SD_RC[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F2_SD_RLOW[m](i, j, k, TETRAMER_DIM_A - 1) = F2_SD_RLOW[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F2_SD_RHIGH[m](i, j, k, TETRAMER_DIM_A - 1) = F2_SD_RHIGH[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);

                        F2_SD_R0[m](TETRAMER_DIM_A - 1, j, k, i) = F2_SD_R0[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F2_SD_RC[m](TETRAMER_DIM_A - 1, j, k, i) = F2_SD_RC[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F2_SD_RLOW[m](TETRAMER_DIM_A - 1, j, k, i) = F2_SD_RLOW[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F2_SD_RHIGH[m](TETRAMER_DIM_A - 1, j, k, i) = F2_SD_RHIGH[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                    }

                    // f3
                    for(int m = 0; m < 7; m++) {
                        _excl_s[m](i, j, k, TETRAMER_DIM_A - 1) = _excl_s[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        _excl_r[m](i, j, k, TETRAMER_DIM_A - 1) = _excl_r[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);

                        _excl_s[m](TETRAMER_DIM_A - 1, j, k, i) = _excl_s[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        _excl_r[m](TETRAMER_DIM_A - 1, j, k, i) = _excl_r[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                    }

                    // f4
                    for(int m = 0; m < 21; m++) {
                        // these ifs skip the parameters we haven't changed; saves a few operations, but it is a bit dangerous when developing the model.
                        // if(m == 2 || m == 3 || m == 4 || m == 5) continue;
                        // if(m == 10 || m == 11 || m == 12) continue;
                        // if(m == 13 || m == 14 || m == 16 || m == 17 || m == 18 || m == 19 || m == 20) continue;
                        F4_SD_THETA_A[m](i, j, k, TETRAMER_DIM_A - 1) = F4_SD_THETA_A[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F4_SD_THETA_T0[m](i, j, k, TETRAMER_DIM_A - 1) = F4_SD_THETA_T0[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F4_SD_THETA_TS[m](i, j, k, TETRAMER_DIM_A - 1) = F4_SD_THETA_TS[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);

                        F4_SD_THETA_A[m](TETRAMER_DIM_A - 1, j, k, i) = F4_SD_THETA_A[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F4_SD_THETA_T0[m](TETRAMER_DIM_A - 1, j, k, i) = F4_SD_THETA_T0[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F4_SD_THETA_TS[m](TETRAMER_DIM_A - 1, j, k, i) = F4_SD_THETA_TS[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                    }

                    // f5
                    for(int m = 0; m < 4; m++) {
                        F5_SD_PHI_A[m](i, j, k, TETRAMER_DIM_A - 1) = F5_SD_PHI_A[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);
                        F5_SD_PHI_XS[m](i, j, k, TETRAMER_DIM_A - 1) = F5_SD_PHI_XS[m].get_average_par(i, j, k, TETRAMER_DIM_A - 1);

                        F5_SD_PHI_A[m](TETRAMER_DIM_A - 1, j, k, i) = F5_SD_PHI_A[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                        F5_SD_PHI_XS[m](TETRAMER_DIM_A - 1, j, k, i) = F5_SD_PHI_XS[m].get_average_par(TETRAMER_DIM_A - 1, j, k ,i);
                    }
                }
            }
        }

        // Set enslaved parameters
        // Enslaved = set by continuity and differentiability
        for(int i = 0; i < TETRAMER_DIM_A; i++) {
            for(int j = 0; j < TETRAMER_DIM_B - 1; j++) {
                for(int k = 0; k < TETRAMER_DIM_B - 1; k++) {
                    for(int l = 0; l < TETRAMER_DIM_A; l++) {
                        // delta2
                        _fene_delta2_SD(i, j, k, l) = SQR(_fene_delta_SD(i, j, k, l));

                        // f1
                        F1_SD_RLOW[0](i, j, k, l) = F1_SD_R0[0](i, j, k, l)-0.06;
                        F1_SD_RHIGH[0](i, j, k, l) = F1_SD_R0[0](i, j, k, l)+0.3;
                        F1_SD_RC[0](i, j, k, l) = F1_SD_R0[0](i, j, k, l)+0.35;

                        F1_SD_RLOW[1](i, j, k, l) = F1_SD_R0[1](i, j, k, l)-0.08;
                        F1_SD_RHIGH[1](i, j, k, l) = F1_SD_R0[1](i, j, k, l)+0.35;
                        F1_SD_RC[1](i, j, k, l) = F1_SD_R0[1](i, j, k, l)+0.5;

                        number term1, term2, term3;
                        for(int m = 0; m < 2; m++) {

                            term1 = exp(-F1_SD_A[m](i, j, k, l) * (F1_SD_RLOW[m](i, j, k, l) - F1_SD_R0[m](i, j, k, l)));
                            term2 = exp(-F1_SD_A[m](i, j, k, l) * (F1_SD_RC[m](i, j, k, l) - F1_SD_R0[m](i, j, k, l)));
                            term3 = exp(-F1_SD_A[m](i, j, k, l) * (F1_SD_RHIGH[m](i, j, k, l) - F1_SD_R0[m](i, j, k, l)));

                            F1_SD_BLOW[m](i, j, k, l) = SQR(F1_SD_A[m](i, j, k, l) * term1 * (1 - term1)) / (SQR(1 - term1) - SQR(1 - term2));
                            F1_SD_BHIGH[m](i, j, k, l) = SQR(F1_SD_A[m](i, j, k, l) * term3 * (1 - term3)) / (SQR(1 - term3) - SQR(1 - term2));

                            F1_SD_RCLOW[m](i, j, k, l) = F1_SD_RLOW[m](i, j, k, l) - F1_SD_A[m](i, j, k, l) / F1_SD_BLOW[m](i, j, k, l) * (term1 * (1 - term1));
                            F1_SD_RCHIGH[m](i, j, k, l) = F1_SD_RHIGH[m](i, j, k, l) - F1_SD_A[m](i, j, k, l) / F1_SD_BHIGH[m](i, j, k, l) * (term3 * (1 - term3));

                            F1_SD_SHIFT[m](i, j, k, l) = F1_SD_EPS[m](i, j, k, l) * SQR(1 - exp(-(F1_SD_RC[m](i, j, k, l) - F1_SD_R0[m](i, j, k, l)) * F1_SD_A[m](i, j, k, l)));
                        }
                        // f2
                        for(int m = 0; m < 4; m++) {

                            F2_SD_RLOW[m](i, j, k, l) = F2_SD_R0[m](i, j, k, l)-0.08;
                            F2_SD_RHIGH[m](i, j, k, l) = F2_SD_R0[m](i, j, k, l)+0.08;
                            F2_SD_RC[m](i, j, k, l) = F2_SD_R0[m](i, j, k, l)+0.1;

                            term1 = F2_SD_RLOW[m](i, j, k, l) - F2_SD_R0[m](i, j, k, l);
                            term2 = F2_SD_RHIGH[m](i, j, k, l) - F2_SD_R0[m](i, j, k, l);
                            term3 = F2_SD_RC[m](i, j, k, l) - F2_SD_R0[m](i, j, k, l);

                            F2_SD_RCLOW[m](i, j, k, l) = F2_SD_RLOW[m](i, j, k, l) - term1 + SQR(term3) / term1;
                            F2_SD_RCHIGH[m](i, j, k, l) = F2_SD_RHIGH[m](i, j, k, l) - term2 + SQR(term3) / term2;

                            F2_SD_BLOW[m](i, j, k, l) = -0.5 * term1 / (F2_SD_RCLOW[m](i, j, k, l) - F2_SD_RLOW[m](i, j, k, l));
                            F2_SD_BHIGH[m](i, j, k, l) = -0.5 * term2 / (F2_SD_RCHIGH[m](i, j, k, l) - F2_SD_RHIGH[m](i, j, k, l));
                        }

                        // f3
                        for(int m = 0; m < 7; m++) {
                            number tmp = SQR(_excl_s[m](i, j, k, l) / _excl_r[m](i, j, k, l));
                            term1 = tmp * tmp * tmp;
                            term2 = 4. * (SQR(term1) - term1);
                            term3 = 12. / _excl_r[m](i, j, k, l) * (2. * SQR(term1) - term1);

                            _excl_rc[m](i, j, k, l) = term2 / term3 + _excl_r[m](i, j, k, l);
                            _excl_b[m](i, j, k, l) = SQR(term3) / term2;
                        }

                        // f4
                        for(int m = 0; m < 21; m++) {
                            F4_SD_THETA_TS[m](i, j, k, l) = sqrt(0.81225/F4_SD_THETA_A[m](i, j, k, l));
                            F4_SD_THETA_TC[m](i, j, k, l) = 1. / F4_SD_THETA_A[m](i, j, k, l) / F4_SD_THETA_TS[m](i, j, k, l);
                            F4_SD_THETA_B[m](i, j, k, l) = F4_SD_THETA_A[m](i, j, k, l) * F4_SD_THETA_TS[m](i, j, k, l) / (F4_SD_THETA_TC[m](i, j, k, l) - F4_SD_THETA_TS[m](i, j, k, l));
                        }
                        // f5
                        for(int m = 0; m < 4; m++) {
                            term1 = 1. - F5_SD_PHI_A[m](i, j, k, l) * SQR(F5_SD_PHI_XS[m](i, j, k, l));
                            term2 = F5_SD_PHI_A[m](i, j, k, l) * F5_SD_PHI_XS[m](i, j, k, l);

                            F5_SD_PHI_XC[m](i, j, k, l) = term1 / term2 + F5_SD_PHI_XS[m](i, j, k, l);
                            F5_SD_PHI_B[m](i, j, k, l) = SQR(term2) / term1;
                        }
                    }
                }
            }
        }

        // updating _mbf_fmax after accounting for new value of fene_delta. In general this is only applied if max_backbone_force is given

        if(_use_mbf) {
            for(int i = 0; i < TETRAMER_DIM_A; i++) {
                for(int j = 0; j < TETRAMER_DIM_B; j++) {
                    for(int k = 0; k < TETRAMER_DIM_B; k++) {
                        for(int l = 0; l < TETRAMER_DIM_A; l++) {
                            _mbf_xmax_SD(i, j, k, l) = (-_fene_eps + sqrt(_fene_eps * _fene_eps + 4.f * _mbf_fmax * _mbf_fmax * _fene_delta2_SD(i, j, k, l))) / (2.f * _mbf_fmax);
                            // if we use mbf, we should tell the user
                            OX_LOG(Logger::LOG_INFO, "Overwriting mbf_xmax %c %c %c %c to %g", Utils::encode_base(i), Utils::encode_base(k), Utils::encode_base(l), Utils::encode_base(j), _mbf_xmax_SD(i, j, k, l));
                        }
                    }
                }
            }
        }

        // We don't have to change this, since we are not touching the unbonded excluded volume.
        // Fixing rcut (cut distance for nonbonded interaction computation)
        // We loop through all bases combinations, and pick the largest cutoff distance possible.
        // Can be made more efficient if we define an SD version and use the right rcut depending on the types of interacting particles,
        // But might not work in VMMC. (Andrea: don't know, might have to check. Also, how much would we gain?).
        number rcutbase;
        number rcutback;
        number rcutbase_max = 0.;
        number rcutback_max = 0.;
        // QUESTION: these three arrays are identical... can we avoid the duplication?
        number pb1[4] = {POS_MM_BACK1_A, POS_MM_BACK1_G, POS_MM_BACK1_C, POS_MM_BACK1_T};
        number pb2[4] = {POS_MM_BACK1_A, POS_MM_BACK1_G, POS_MM_BACK1_C, POS_MM_BACK1_T};
        number pba[4] = {POS_MM_BACK1_A, POS_MM_BACK1_G, POS_MM_BACK1_C, POS_MM_BACK1_T};
        for(int i = 0; i < 6; i++) {
            for(int j = 0; j < 4; j++) {
                for(int k = 0; k < 4; k++) {
                    for(int l = 0; l < 6; l++) {
                        for(int m = 0; m < 4; m++) {
                            for(int n = 0; n < 4; n++) {
                                rcutback = sqrt((pb1[m]) * (pb1[m]) + (pb2[m]) * (pb2[m])) + sqrt((pb1[n]) * (pb1[n]) + (pb2[n]) * (pb2[n])) + _excl_rc[0](i, j, k, l);
                                if(rcutback > rcutback_max) rcutback_max = rcutback;
                                rcutbase = fabs(pba[n]) + fabs(pba[m]) + F1_SD_RCHIGH[0](i, j, k, l);
                                if(rcutbase > rcutbase_max) rcutbase_max = rcutbase;
                            }
                        }
                    }
                }
            }
        }

        _rcut = fmax(rcutback_max, rcutbase_max);
        _sqr_rcut = SQR(_rcut);

        number lambda = _debye_huckel_lambdafactor * sqrt(_T / 0.1) / sqrt(_salt_concentration);
        _minus_kappa = -1.0 / lambda;

        // these are just for convenience for the smoothing parameter computation
        number exs = _debye_huckel_RHIGH;
        number q = _debye_huckel_prefactor;
        number la = lambda;

        // compute the some smoothing parameters
        _debye_huckel_B = -(exp(-exs / la) * q * q * (exs + la) * (exs + la)) / (-4. * exs * exs * exs * la * la * q);
        _debye_huckel_RC = exs * (q * exs + 3. * q * la) / (q * (exs + la));

        // note: in oxdna3 POS_MM_BACK1 and POS_MM_BACK2 are the same for all nucleotides
        number debyecut = 2.0 * sqrt(SQR(POS_MM_BACK1) + SQR(POS_MM_BACK2)) + _debye_huckel_RC;
        // the cutoff radius for the potential should be the larger of rcut and debyecut
        if(debyecut > _rcut) {
            _rcut = debyecut;
            _sqr_rcut = SQR(_rcut);
        }

        OX_LOG(Logger::LOG_DEBUG,"Debye-Huckel parameters: Q=%f, lambda_0=%f, lambda=%f, r_high=%f, cutoff=%f", _debye_huckel_prefactor, _debye_huckel_lambdafactor, lambda, _debye_huckel_RHIGH, _rcut);
        OX_LOG(Logger::LOG_DEBUG,"Debye-Huckel parameters: debye_huckel_RC=%e, debye_huckel_B=%e", _debye_huckel_RC, _debye_huckel_B);
        OX_LOG(Logger::LOG_INFO,"The Debye length at this temperature (%lf) and salt concentration (%lf) is %f", _T, _salt_concentration, lambda);

        // build mesh for the f4s
        for(int int_type = 0; int_type < 21; int_type++) {
            // the order of the interpolation interval extremes is reversed,
            // due to the cosine being monotonically decreasing with increasing x
            for(int i = 0; i < TETRAMER_DIM_A; i++) {
                for(int j = 0; j < TETRAMER_DIM_B; j++) {
                    for(int k = 0; k < TETRAMER_DIM_B; k++) {
                        for(int l = 0; l < TETRAMER_DIM_A; l++) {
                            int points = MESH_F4_SD_POINTS[int_type];
                            number upplimit = cos(fmax(0, F4_SD_THETA_T0[int_type](i, j, k, l) - F4_SD_THETA_TC[int_type](i, j, k, l)));
                            number lowlimit = cos(fmin(PI, F4_SD_THETA_T0[int_type](i, j, k, l) + F4_SD_THETA_TC[int_type](i, j, k, l)));

                            int tmp_args[5] = {int_type, i, j, k, l};
                            _mesh_f4_SD[int_type][i][j][k][l].build([this](number x, void *args) { return this->_fakef4_SD(x, args); }, [this](number x, void *args) { return _fakef4D_SD(x, args); }, (void *)tmp_args, points, lowlimit, upplimit);
                            assert(lowlimit < upplimit);
                        }
                    }
                }
            }
        }
    }
}

number DNA3Interaction::_bonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(!_check_bonded_neighbour(&p, &q, compute_r)) {
        return (number)0.f;
    }

    if(compute_r) {
        _computed_r = q->pos - p->pos;
    }

    // BASE-BASE
    LR_vector force;
    LR_vector torqueq(0, 0, 0);
    LR_vector torquep(0, 0, 0);

    int type_n3_2 = NO_TYPE;
    int type_n5_2 = NO_TYPE;

    if(q->n3 != P_VIRTUAL) type_n3_2 = q->n3->type;
    if(p->n5 != P_VIRTUAL) type_n5_2 = p->n5->type;

    LR_vector rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
    number s = _excl_s[4](type_n3_2, q->type, p->type, type_n5_2);
    number rstar = _excl_r[4](type_n3_2, q->type, p->type, type_n5_2);
    number b = _excl_b[4](type_n3_2, q->type, p->type, type_n5_2);
    number rc = _excl_rc[4](type_n3_2, q->type, p->type, type_n5_2);
    number energy = _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
        torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);
    }

    // P-BASE vs. Q-BACK
    rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BASE];
    s = _excl_s[5](type_n3_2, q->type, p->type, type_n5_2);
    rstar = _excl_r[5](type_n3_2, q->type, p->type, type_n5_2);
    b = _excl_b[5](type_n3_2, q->type, p->type, type_n5_2);
    rc = _excl_rc[5](type_n3_2, q->type, p->type, type_n5_2);
    energy += _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
        torqueq += q->int_centers[DNANucleotide::BACK].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);
    }

    // P-BACK vs. Q-BASE
    rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BACK];
    s = _excl_s[6](type_n3_2, q->type, p->type, type_n5_2);
    rstar = _excl_r[6](type_n3_2, q->type, p->type, type_n5_2);
    b = _excl_b[6](type_n3_2, q->type, p->type, type_n5_2);
    rc = _excl_rc[6](type_n3_2, q->type, p->type, type_n5_2);
    energy += _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep -= p->int_centers[DNANucleotide::BACK].cross(force);
        torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);

        // we need torques in the reference system of the particle
        p->torque += p->orientationT * torquep;
        q->torque += q->orientationT * torqueq;
    }

    return energy;
}

number DNA3Interaction::_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->is_bonded(q)) {
        return (number)0.f;
    }

    LR_vector force(0, 0, 0);
    LR_vector torquep(0, 0, 0);
    LR_vector torqueq(0, 0, 0);

    // for non bonded we don't look at the flanking; we just take a length that can accomodate all crst and hydr
    // otherwise we would have to introduce a ton of parameters
    int type_n3_2 = NO_TYPE;
    int type_n5_2 = NO_TYPE;

    // BASE-BASE
    LR_vector rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
    number s = _excl_s[1](type_n3_2, q->type, p->type, type_n5_2);
    number rstar = _excl_r[1](type_n3_2, q->type, p->type, type_n5_2);
    number b = _excl_b[1](type_n3_2, q->type, p->type, type_n5_2);
    number rc = _excl_rc[1](type_n3_2, q->type, p->type, type_n5_2);
    number energy = _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep = -p->int_centers[DNANucleotide::BASE].cross(force);
        torqueq = q->int_centers[DNANucleotide::BASE].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);
    }

    // P-BACK vs. Q-BASE
    rcenter = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BACK];
    s = _excl_s[3](type_n3_2, q->type, p->type, type_n5_2);
    rstar = _excl_r[3](type_n3_2, q->type, p->type, type_n5_2);
    b = _excl_b[3](type_n3_2, q->type, p->type, type_n5_2);
    rc = _excl_rc[3](type_n3_2, q->type, p->type, type_n5_2);
    energy += _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep += -p->int_centers[DNANucleotide::BACK].cross(force);
        torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);
    }

    // P-BASE vs. Q-BACK
    rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BASE];
    s = _excl_s[2](type_n3_2, q->type, p->type, type_n5_2);
    rstar = _excl_r[2](type_n3_2, q->type, p->type, type_n5_2);
    b = _excl_b[2](type_n3_2, q->type, p->type, type_n5_2);
    rc = _excl_rc[2](type_n3_2, q->type, p->type, type_n5_2);
    energy += _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep += -p->int_centers[DNANucleotide::BASE].cross(force);
        torqueq += q->int_centers[DNANucleotide::BACK].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);
    }

    // BACK-BACK
    rcenter = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
    s = _excl_s[0](type_n3_2, q->type, p->type, type_n5_2);
    rstar = _excl_r[0](type_n3_2, q->type, p->type, type_n5_2);
    b = _excl_b[0](type_n3_2, q->type, p->type, type_n5_2);
    rc = _excl_rc[0](type_n3_2, q->type, p->type, type_n5_2);
    energy += _repulsive_lj(rcenter, force, s, rstar, b, rc, update_forces);

    if(update_forces) {
        torquep += -p->int_centers[DNANucleotide::BACK].cross(force);
        torqueq += q->int_centers[DNANucleotide::BACK].cross(force);

        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);

        // we need torques in the reference system of the particle
        p->torque += p->orientationT * torquep;
        q->torque += q->orientationT * torqueq;
    }

    return energy;
}

number DNA3Interaction::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(!_check_bonded_neighbour(&p, &q, compute_r)) {
        return (number)0.f;
    }

    if(compute_r) {
        _computed_r = q->pos - p->pos;
    }

    int type_n3_2 = (q->n3 != P_VIRTUAL) ? q->n3->type : NO_TYPE;
    int type_n5_2 = (p->n5 != P_VIRTUAL) ? p->n5->type : NO_TYPE;

    // These two conditions are needed if we impose complementarity symmetry in bonded interaction.
    // if(q->n3 == P_VIRTUAL && p->n5 != P_VIRTUAL) type_n3_2 = type_n5_2;
    // if(p->n5 == P_VIRTUAL && q->n3 != P_VIRTUAL) type_n5_2 = type_n3_2;

    // std::cout << "ptypes: " << q->get_index() << "," << p->get_index() << ", " << type_n3_2 << " " << q->type << " " << p->type << " " <<  type_n5_2 << std::endl;

    LR_vector rback = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
    number rbackmod = rback.module();
    number rbackr0 = rbackmod - _fene_r0_SD(type_n3_2, q->type, p->type, type_n5_2);

    LR_vector force;
    number energy;

    number mbf_xmax = _mbf_xmax_SD(type_n3_2, q->type, p->type, type_n5_2);
    number fene_delta2 = _fene_delta2_SD(type_n3_2, q->type, p->type, type_n5_2);
    if(_use_mbf && fabs(rbackr0) > mbf_xmax) {
        // we use the "log"  potential
        number fene_xmax = -(_fene_eps / 2.f) * log(1.f - SQR(mbf_xmax) /fene_delta2);
        number long_xmax = (_mbf_fmax - _mbf_finf) * mbf_xmax * log(mbf_xmax) + _mbf_finf * mbf_xmax;
        energy = (_mbf_fmax - _mbf_finf) * mbf_xmax * log(fabs(rbackr0)) + _mbf_finf * fabs(rbackr0) - long_xmax + fene_xmax;
        if(update_forces)
            force = rback * (-copysign(1.f, rbackr0) * ((_mbf_fmax - _mbf_finf) * mbf_xmax / fabs(rbackr0) + _mbf_finf) / rbackmod);
    }
    else {
        // we check whether we ended up OUTSIDE of the FENE range
        if(fabs(rbackr0) > _fene_delta_SD(type_n3_2, q->type, p->type, type_n5_2) - DBL_EPSILON) {
            std::cout << "Fene too large: " << q->get_index() << "," << p->get_index() << " type: " << type_n3_2 << q->type << p->type << type_n5_2 << " " << rbackmod << " " << _fene_r0_SD(type_n3_2, q->type, p->type, type_n5_2) << " " << _fene_delta_SD(type_n3_2, q->type, p->type, type_n5_2) << " " << fene_delta2 << std::endl;
            if(update_forces && !_allow_broken_fene) {
                throw oxDNAException("(DNAInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf), %lf, %lf, %lf", p->index, q->index, fabs(rbackr0), fabs(rbackmod), fabs(_fene_r0_SD(type_n3_2, q->type, p->type, type_n5_2)), fabs(_fene_delta_SD(type_n3_2, q->type, p->type, type_n5_2)));
            }
            return (number)(1.e12);
        }

        // if not, we just do the right thing
        energy = -(_fene_eps / 2.f) * log(1.f - SQR(rbackr0) / fene_delta2);
        if(update_forces) {
            force = rback * (-(_fene_eps * rbackr0 / (fene_delta2 - SQR(rbackr0))) / rbackmod);
        }
    }

    if(update_forces) {
        p->force -= force;
        q->force += force;

        _update_stress_tensor(_computed_r, force);

        p->torque -= p->orientationT * p->int_centers[DNANucleotide::BACK].cross(force);
        q->torque += q->orientationT * q->int_centers[DNANucleotide::BACK].cross(force);
    }

    return energy;
}

number DNA3Interaction::_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(!_check_bonded_neighbour(&p, &q, compute_r)) {
        return (number)0.f;
    }

    LR_vector &a1 = p->orientationT.v1;
    LR_vector &a2 = p->orientationT.v2;
    LR_vector &a3 = p->orientationT.v3;
    LR_vector &b1 = q->orientationT.v1;
    LR_vector &b2 = q->orientationT.v2;
    LR_vector &b3 = q->orientationT.v3;

    if(compute_r) {
        _computed_r = q->pos - p->pos;
    }

    // This is the position the backbone would have with major-minor grooves the same width.
    // We need to do this to implement different major-minor groove widths because rback is
    // used as a reference point for things that have nothing to do with the actual backbone
    // position (in this case, the stacking interaction).
    LR_vector rbackref = _computed_r + b1 * POS_BACK - a1 * POS_BACK;
    number rbackrefmod = rbackref.module();

    LR_vector rstack = _computed_r + q->int_centers[DNANucleotide::STACK] - p->int_centers[DNANucleotide::STACK];
    number rstackmod = rstack.module();
    LR_vector rstackdir = rstack / rstackmod;

    number cost4 = a3 * b3;
    number cost5 = a3 * rstackdir;
    number cost6 = -b3 * rstackdir;
    number cosphi1 = a2 * rbackref / rbackrefmod;
    number cosphi2 = b2 * rbackref / rbackrefmod;

    // functions and their derivatives needed for energies and forces
    int type_n3_2 = NO_TYPE;
    int type_n5_2 = NO_TYPE;

    if(q->n3 != P_VIRTUAL) type_n3_2 = q->n3->type;
    if(p->n5 != P_VIRTUAL) type_n5_2 = p->n5->type;

    // These two conditions are needed if we impose complementarity symmetry in bonded interaction.
    // if(q->n3 == P_VIRTUAL && p->n5 != P_VIRTUAL) type_n3_2 = type_n5_2;
    // if(p->n5 == P_VIRTUAL && q->n3 != P_VIRTUAL) type_n5_2 = type_n3_2;

    number f1 = _f1_SD(rstackmod, STCK_F1, type_n3_2, q->type, p->type, type_n5_2);
    number f4t4 = _custom_f4_SD(cost4, STCK_F4_THETA4, type_n3_2, q->type, p->type, type_n5_2);
    number f4t5 = _custom_f4_SD(-cost5, STCK_F4_THETA5, type_n3_2, q->type, p->type, type_n5_2);
    number f4t6 = _custom_f4_SD(cost6, STCK_F4_THETA6, type_n3_2, q->type, p->type, type_n5_2);
    number f5phi1 = _f5_SD(cosphi1, STCK_F5_PHI1, type_n3_2, q->type, p->type, type_n5_2);
    number f5phi2 = _f5_SD(cosphi2, STCK_F5_PHI2, type_n3_2, q->type, p->type, type_n5_2);

    number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

    if(update_forces && energy != (number)0.f) {
        LR_vector torquep(0, 0, 0);
        LR_vector torqueq(0, 0, 0);

        // these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
        // we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
        number f1D = _f1D_SD(rstackmod, STCK_F1, type_n3_2, q->type, p->type, type_n5_2);
        number f4t4Dsin = -_custom_f4D_SD(cost4, STCK_F4_THETA4, type_n3_2, q->type, p->type, type_n5_2);
        number f4t5Dsin = -_custom_f4D_SD(-cost5, STCK_F4_THETA5, type_n3_2, q->type, p->type, type_n5_2);
        number f4t6Dsin = -_custom_f4D_SD(cost6, STCK_F4_THETA6, type_n3_2, q->type, p->type, type_n5_2);
        number f5phi1D = _f5D_SD(cosphi1, STCK_F5_PHI1, type_n3_2, q->type, p->type, type_n5_2);
        number f5phi2D = _f5D_SD(cosphi2, STCK_F5_PHI2, type_n3_2, q->type, p->type, type_n5_2);

        // RADIAL
        LR_vector force = -rstackdir * (f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2);

        // THETA 5
        force += -(a3 - rstackdir * cost5) * (f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 / rstackmod);

        // THETA 6
        force += -(b3 + rstackdir * cost6) * (f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 / rstackmod);

        // COS PHI 1
        // here particle p is referred to using the a while particle q is referred with the b
        number gamma = POS_STACK - POS_BACK;
        number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;

        number ra2 = rstackdir * a2;
        number ra1 = rstackdir * a1;
        number rb1 = rstackdir * b1;
        number a2b1 = a2 * b1;

        number parentesi = rstackmod * ra2 - a2b1 * gamma;

        number dcosphi1dr = (SQR(rstackmod) * ra2 - ra2 * SQR(rbackrefmod) - rstackmod * (a2b1 + ra2 * (-ra1 + rb1)) * gamma + a2b1 * (-ra1 + rb1) * SQR(gamma)) / rbackrefmodcub;
        number dcosphi1dra1 = rstackmod * gamma * parentesi / rbackrefmodcub;
        number dcosphi1dra2 = -rstackmod / rbackrefmod;
        number dcosphi1drb1 = -rstackmod * gamma * parentesi / rbackrefmodcub;
        number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackrefmodcub;
        number dcosphi1da2b1 = gamma / rbackrefmod;

        // this force part has a minus because of the definition of cos(phi1)
        // which is not consistent with the derivatives above (i.e. all the
        // derivatives should have a minus sign in front and the force below shouldn't)
        number force_part_phi1 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1D * f5phi2;

        force += -(rstackdir * dcosphi1dr +
                   ((a2 - rstackdir * ra2) * dcosphi1dra2 +
                    (a1 - rstackdir * ra1) * dcosphi1dra1 +
                    (b1 - rstackdir * rb1) * dcosphi1drb1) /
                       rstackmod) *
                 force_part_phi1;

        // COS PHI 2
        // here particle p -> b, particle q -> a
        ra2 = rstackdir * b2;
        ra1 = rstackdir * b1;
        rb1 = rstackdir * a1;
        a2b1 = b2 * a1;

        parentesi = rstackmod * ra2 + a2b1 * gamma;

        number dcosphi2dr = (parentesi * (rstackmod + (rb1 - ra1) * gamma) - ra2 * SQR(rbackrefmod)) / rbackrefmodcub;
        number dcosphi2dra1 = -rstackmod * gamma * (rstackmod * ra2 + a2b1 * gamma) / rbackrefmodcub;
        number dcosphi2dra2 = -rstackmod / rbackrefmod;
        number dcosphi2drb1 = rstackmod * gamma * parentesi / rbackrefmodcub;
        number dcosphi2da1b1 = -SQR(gamma) * parentesi / rbackrefmodcub;
        number dcosphi2da2b1 = -gamma / rbackrefmod;

        // this force part has a minus because of the definition of cos(phi2) which is not consistent with the derivatives
        // above (i.e. all the derivatives should have a minus sign in front and the force below shouldn't)
        number force_part_phi2 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2D;

        force += -force_part_phi2 * (rstackdir * dcosphi2dr +
                                     ((b2 - rstackdir * ra2) * dcosphi2dra2 +
                                      (b1 - rstackdir * ra1) * dcosphi2dra1 +
                                      (a1 - rstackdir * rb1) * dcosphi2drb1) /
                                         rstackmod);

        // Add the contribution from all the forces to the stored particles' forces
        p->force -= force;
        q->force += force;

        torquep -= p->int_centers[DNANucleotide::STACK].cross(force);
        torqueq += q->int_centers[DNANucleotide::STACK].cross(force);

        // handle the part on theta4
        LR_vector t4dir = b3.cross(a3);
        number torquemod = f1 * f4t4Dsin * f4t5 * f4t6 * f5phi1 * f5phi2;

        torquep -= t4dir * torquemod;
        torqueq += t4dir * torquemod;

        // handle the part on theta5
        LR_vector t5dir = rstackdir.cross(a3);
        torquemod = -f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;

        torquep -= t5dir * torquemod;

        // handle the part on theta6
        LR_vector t6dir = rstackdir.cross(b3);
        torquemod = f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;

        torqueq += t6dir * torquemod;

        // PHI 1
        torquep += rstackdir.cross(a2) * force_part_phi1 * dcosphi1dra2 +
                   rstackdir.cross(a1) * force_part_phi1 * dcosphi1dra1;
        torqueq += rstackdir.cross(b1) * force_part_phi1 * dcosphi1drb1;

        LR_vector puretorque = a2.cross(b1) * force_part_phi1 * dcosphi1da2b1 +
                               a1.cross(b1) * force_part_phi1 * dcosphi1da1b1;

        torquep -= puretorque;
        torqueq += puretorque;

        // PHI 2
        torquep += rstackdir.cross(a1) * force_part_phi2 * dcosphi2drb1;
        torqueq += rstackdir.cross(b2) * force_part_phi2 * dcosphi2dra2 +
                   rstackdir.cross(b1) * force_part_phi2 * dcosphi2dra1;

        puretorque = a1.cross(b2) * force_part_phi2 * dcosphi2da2b1 +
                     a1.cross(b1) * force_part_phi2 * dcosphi2da1b1;

        torquep -= puretorque;
        torqueq += puretorque;

        // we need torques in the reference system of the particle
        p->torque += p->orientationT * torquep;
        q->torque += q->orientationT * torqueq;
    }

    return energy;
}

number DNA3Interaction::_hydrogen_bonding(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->is_bonded(q)) {
        return (number)0.f;
    }

    // true if p and q are Watson-Crick-like pairs
    bool is_pair = (q->btype + p->btype == 3);
    number hb_multi = (abs(q->btype) >= 300 && abs(p->btype) >= 300) ? _hb_multiplier : 1.f;

    LR_vector rhydro = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
    number rhydromod = rhydro.module();
    number energy = (number)0.f;
    if(is_pair && F1_SD_RCLOW[HYDR_F1](0, q->type, p->type, 0) < rhydromod && rhydromod < F1_SD_RCHIGH[HYDR_F1](0, q->type, p->type, 0)) {
        // vector, versor and magnitude of the base-base separation
        LR_vector rhydrodir = rhydro / rhydromod;

        // particle axes according to Allen's paper
        LR_vector &a1 = p->orientationT.v1;
        LR_vector &a3 = p->orientationT.v3;
        LR_vector &b1 = q->orientationT.v1;
        LR_vector &b3 = q->orientationT.v3;

        // angles involved in the HB interaction
        number cost1 = -a1 * b1;
        number cost2 = -b1 * rhydrodir;
        number cost3 = a1 * rhydrodir;

        number cost4 = a3 * b3;
        number cost7 = -b3 * rhydrodir;
        number cost8 = a3 * rhydrodir;

        // functions called at their relevant arguments
        number f1 = hb_multi * _f1_SD(rhydromod, HYDR_F1, 0, q->type, p->type, 0);
        number f4t1 = _custom_f4_SD(cost1, HYDR_F4_THETA1, 0, q->type, p->type, 0);
        number f4t2 = _custom_f4_SD(cost2, HYDR_F4_THETA2, 0, q->type, p->type, 0);
        number f4t3 = _custom_f4_SD(cost3, HYDR_F4_THETA3, 0, q->type, p->type, 0);

        number f4t4 = _custom_f4_SD(cost4, HYDR_F4_THETA4, 0, q->type, p->type, 0);
        number f4t7 = _custom_f4_SD(cost7, HYDR_F4_THETA7, 0, q->type, p->type, 0);
        number f4t8 = _custom_f4_SD(cost8, HYDR_F4_THETA8, 0, q->type, p->type, 0);

        energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

        //  makes sense, since the above functions may return 0. exactly
        if(update_forces && energy != 0.) {
            LR_vector force(0, 0, 0);
            LR_vector torquep(0, 0, 0);
            LR_vector torqueq(0, 0, 0);

            // derivatives called at the relevant arguments
            number f1D = hb_multi * _f1D_SD(rhydromod, HYDR_F1, 0, q->type, p->type, 0);
            number f4t1Dsin = _custom_f4D_SD(cost1, HYDR_F4_THETA1, 0, q->type, p->type, 0);
            number f4t2Dsin = _custom_f4D_SD(cost2, HYDR_F4_THETA2, 0, q->type, p->type, 0);
            number f4t3Dsin = -_custom_f4D_SD(cost3, HYDR_F4_THETA3, 0, q->type, p->type, 0);

            number f4t4Dsin = -_custom_f4D_SD(cost4, HYDR_F4_THETA4, 0, q->type, p->type, 0);
            number f4t7Dsin = _custom_f4D_SD(cost7, HYDR_F4_THETA7, 0, q->type, p->type, 0);
            number f4t8Dsin = -_custom_f4D_SD(cost8, HYDR_F4_THETA8, 0, q->type, p->type, 0);

            // RADIAL PART
            force = -rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

            // THETA4; t4 = LRACOS (a3 * b3);
            LR_vector dir = a3.cross(b3);
            number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

            torquep -= dir * torquemod;
            torqueq += dir * torquemod;

            // THETA1; t1 = LRACOS (-a1 * b1);
            dir = a1.cross(b1);
            torquemod = -f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

            torquep -= dir * torquemod;
            torqueq += dir * torquemod;

            // THETA2; t2 = LRACOS (-b1 * rhydrodir);
            number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
            force += (b1 + rhydrodir * cost2) * (fact / rhydromod);
            dir = rhydrodir.cross(b1);

            torqueq += -dir * fact;

            // THETA3; t3 = LRACOS (a1 * rhydrodir);
            force += (a1 - rhydrodir * cost3) * (f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);

            LR_vector t3dir = rhydrodir.cross(a1);
            torquemod = -f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

            torquep += t3dir * torquemod;

            // THETA7; t7 = LRACOS (-rhydrodir * b3);
            force += (b3 + rhydrodir * cost7) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rhydromod);

            LR_vector t7dir = rhydrodir.cross(b3);
            torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

            torqueq += t7dir * torquemod;

            // THETA 8; t8 = LRACOS (rhydrodir * a3);
            force += (a3 - rhydrodir * cost8) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);

            LR_vector t8dir = rhydrodir.cross(a3);
            torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

            torquep += t8dir * torquemod;

            // final update of forces and torques for HB
            p->force -= force;
            q->force += force;

            torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
            torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

            // we need torques in the reference system of the particle
            p->torque += p->orientationT * torquep;
            q->torque += q->orientationT * torqueq;
        }
    }

    return energy;
}

number DNA3Interaction::_cross_stacking(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    if(p->is_bonded(q)) {
        return (number)0.f;
    }

    LR_vector rcstack = _computed_r + q->int_centers[DNANucleotide::BASE] - p->int_centers[DNANucleotide::BASE];
    number rcstackmod = rcstack.module();
    number energy = (number)0.f;

    LR_vector rcstackdir = rcstack / rcstackmod;

    // particle axes according to Allen's paper
    LR_vector &a1 = p->orientationT.v1;
    LR_vector &a3 = p->orientationT.v3;
    LR_vector &b1 = q->orientationT.v1;
    LR_vector &b3 = q->orientationT.v3;

    number cost7 = -b3 * rcstackdir;
    number cost8 = a3 * rcstackdir;

    int type_p_n3 = NO_TYPE;
    int type_q_n3 = NO_TYPE;

    if(p->n3 != P_VIRTUAL) type_p_n3 = p->n3->type;
    if(q->n3 != P_VIRTUAL) type_q_n3 = q->n3->type;

    // if(q->n5 == P_VIRTUAL && p->n5 != P_VIRTUAL) type_n3_2_33 = type_n5_2_33;
    // if(p->n5 == P_VIRTUAL && q->n5 != P_VIRTUAL) type_n5_2_33 = type_n3_2_33;

    int type_p_n5 = NO_TYPE;
    int type_q_n5 = NO_TYPE;

    if(p->n5 != P_VIRTUAL) type_p_n5 = p->n5->type;
    if(q->n5 != P_VIRTUAL) type_q_n5 = q->n5->type;

    // if(q->n3 == P_VIRTUAL && p->n3 != P_VIRTUAL) type_n3_2_55 = type_n5_2_55;
    // if(p->n3 == P_VIRTUAL && q->n3 != P_VIRTUAL) type_n5_2_55 = type_n3_2_55;

    if((cost7 > 0 && cost8 > 0 && F2_SD_RCLOW[CRST_F2_33](type_q_n3, q->type, p->type, type_p_n3) < rcstackmod && rcstackmod < F2_SD_RCHIGH[CRST_F2_33](type_q_n3, q->type, p->type, type_p_n3)) || (cost7 < 0 && cost8 < 0 && F2_SD_RCLOW[CRST_F2_55](type_q_n5, q->type, p->type, type_p_n5) < rcstackmod && rcstackmod < F2_SD_RCHIGH[CRST_F2_55](type_q_n5, q->type, p->type, type_p_n5))) {
        // angles involved in the CRST interaction
        number cost1 = -a1 * b1;
        number cost2 = -b1 * rcstackdir;
        number cost3 = a1 * rcstackdir;
        number cost4 = a3 * b3;

        // functions called at their relevant arguments
        // 3'3' diagonal
        number f2_33 = _f2_SD(rcstackmod, CRST_F2_33, type_q_n3, q->type, p->type, type_p_n3);
        number f4t1_33 = _custom_f4_SD(cost1, CRST_F4_THETA1_33, type_q_n3, q->type, p->type, type_p_n3);
        number f4t2_33 = _custom_f4_SD(cost2, CRST_F4_THETA2_33, type_q_n3, q->type, p->type, type_p_n3);
        number f4t3_33 = _custom_f4_SD(cost3, CRST_F4_THETA3_33, type_q_n3, q->type, p->type, type_p_n3);
        number f4t4_33 = _custom_f4_SD(cost4, CRST_F4_THETA4_33, type_q_n3, q->type, p->type, type_p_n3);
        number f4t7_33 = _custom_f4_SD(cost7, CRST_F4_THETA7_33, type_q_n3, q->type, p->type, type_p_n3);
        number f4t8_33 = _custom_f4_SD(cost8, CRST_F4_THETA8_33, type_q_n3, q->type, p->type, type_p_n3);

        // 5'5' diagonal
        number f2_55 = _f2_SD(rcstackmod, CRST_F2_55, type_q_n5, q->type, p->type, type_p_n5);
        number f4t1_55 = _custom_f4_SD(cost1, CRST_F4_THETA1_55, type_q_n5, q->type, p->type, type_p_n5);
        number f4t2_55 = _custom_f4_SD(cost2, CRST_F4_THETA2_55, type_q_n5, q->type, p->type, type_p_n5);
        number f4t3_55 = _custom_f4_SD(cost3, CRST_F4_THETA3_55, type_q_n5, q->type, p->type, type_p_n5);
        number f4t4_55 = _custom_f4_SD(cost4, CRST_F4_THETA4_55, type_q_n5, q->type, p->type, type_p_n5);
        number f4t7_55 = _custom_f4_SD(cost7, CRST_F4_THETA7_55, type_q_n5, q->type, p->type, type_p_n5);
        number f4t8_55 = _custom_f4_SD(cost8, CRST_F4_THETA8_55, type_q_n5, q->type, p->type, type_p_n5);

        energy = f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

        // makes sense since the above functions can return exactly 0
        if(update_forces && energy != (number)0.f) {
            LR_vector force(0, 0, 0);
            LR_vector torquep(0, 0, 0);
            LR_vector torqueq(0, 0, 0);

            // derivatives called at the relevant arguments
            number f2D_33 = _f2D_SD(rcstackmod, CRST_F2_33, type_q_n3, q->type, p->type, type_p_n3);
            // std::cout << f2D_33 << std::endl;
            number f4t1Dsin_33 = _custom_f4D_SD(cost1, CRST_F4_THETA1_33, type_q_n3, q->type, p->type, type_p_n3);
            number f4t2Dsin_33 = _custom_f4D_SD(cost2, CRST_F4_THETA2_33, type_q_n3, q->type, p->type, type_p_n3);
            number f4t3Dsin_33 = -_custom_f4D_SD(cost3, CRST_F4_THETA3_33, type_q_n3, q->type, p->type, type_p_n3);
            number f4t4Dsin_33 = -_custom_f4D_SD(cost4, CRST_F4_THETA4_33, type_q_n3, q->type, p->type, type_p_n3);
            number f4t7Dsin_33 = _custom_f4D_SD(cost7, CRST_F4_THETA7_33, type_q_n3, q->type, p->type, type_p_n3);
            number f4t8Dsin_33 = -_custom_f4D_SD(cost8, CRST_F4_THETA8_33, type_q_n3, q->type, p->type, type_p_n3);

            number f2D_55 = _f2D_SD(rcstackmod, CRST_F2_55, type_q_n5, q->type, p->type, type_p_n5);
            // std::cout << f2D_55 << std::endl;
            number f4t1Dsin_55 = _custom_f4D_SD(cost1, CRST_F4_THETA1_55, type_q_n5, q->type, p->type, type_p_n5);
            number f4t2Dsin_55 = _custom_f4D_SD(cost2, CRST_F4_THETA2_55, type_q_n5, q->type, p->type, type_p_n5);
            number f4t3Dsin_55 = -_custom_f4D_SD(cost3, CRST_F4_THETA3_55, type_q_n5, q->type, p->type, type_p_n5);
            number f4t4Dsin_55 = -_custom_f4D_SD(cost4, CRST_F4_THETA4_55, type_q_n5, q->type, p->type, type_p_n5);
            number f4t7Dsin_55 = _custom_f4D_SD(cost7, CRST_F4_THETA7_55, type_q_n5, q->type, p->type, type_p_n5);
            number f4t8Dsin_55 = -_custom_f4D_SD(cost8, CRST_F4_THETA8_55, type_q_n5, q->type, p->type, type_p_n5);

            // RADIAL PART
            force = -rcstackdir * ((f2D_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33) + (f2D_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55));

            // THETA1; t1 = LRACOS (-a1 * b1);
            LR_vector t1dir = a1.cross(b1);
            number torquemod = -f2_33 * f4t1Dsin_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1Dsin_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

            torquep -= t1dir * torquemod;
            torqueq += t1dir * torquemod;

            // THETA2; t2 = LRACOS (-b1 * rhydrodir);
            force += (b1 + rcstackdir * cost2) * ((f2_33 * f4t1_33 * f4t2Dsin_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2Dsin_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55) / rcstackmod);

            LR_vector t2dir = rcstackdir.cross(b1);
            torquemod = -f2_33 * f4t1_33 * f4t2Dsin_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1_55 * f4t2Dsin_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

            torqueq += t2dir * torquemod;

            // THETA3; t3 = LRACOS (a1 * rhydrodir);
            force += (a1 - rcstackdir * cost3) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3Dsin_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3Dsin_55 * f4t4_55 * f4t7_55 * f4t8_55) / rcstackmod);

            LR_vector t3dir = rcstackdir.cross(a1);
            torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3Dsin_33 * f4t4_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3Dsin_55 * f4t4_55 * f4t7_55 * f4t8_55;

            torquep += t3dir * torquemod;

            // THETA4; t4 = LRACOS (a3 * b3);
            LR_vector t4dir = a3.cross(b3);
            torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4Dsin_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4Dsin_55 * f4t7_55 * f4t8_55;

            torquep -= t4dir * torquemod;
            torqueq += t4dir * torquemod;

            // THETA7; t7 = LRACOS (-rcsrackir * b3);
            force += (b3 + rcstackdir * cost7) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7Dsin_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7Dsin_55 * f4t8_55) / rcstackmod);
            LR_vector t7dir = rcstackdir.cross(b3);
            torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7Dsin_33 * f4t8_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7Dsin_55 * f4t8_55;

            torqueq += t7dir * torquemod;

            // THETA 8; t8 = LRACOS (rhydrodir * a3);
            force += (a3 - rcstackdir * cost8) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8Dsin_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8Dsin_55) / rcstackmod);

            LR_vector t8dir = rcstackdir.cross(a3);
            torquemod = -f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8Dsin_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8Dsin_55;

            torquep += t8dir * torquemod;

            // final update of forces and torques for CRST
            p->force -= force;
            q->force += force;

            torquep -= p->int_centers[DNANucleotide::BASE].cross(force);
            torqueq += q->int_centers[DNANucleotide::BASE].cross(force);

            // we need torques in the reference system of the particle
            p->torque += p->orientationT * torquep;
            q->torque += q->orientationT * torqueq;
        }
    }

    return energy;
}

number DNA3Interaction::_f1_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;
    if(r < F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
        if(r > F1_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
            val = F1_SD_EPS[type](n3_2, n3_1, n5_1, n5_2) * F1_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
        } else if(r > F1_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            number tmp = 1 - exp(-(r - F1_SD_R0[type](n3_2, n3_1, n5_1, n5_2)) * F1_SD_A[type](n3_2, n3_1, n5_1, n5_2));
            val = F1_SD_EPS[type](n3_2, n3_1, n5_1, n5_2) * SQR(tmp) - F1_SD_SHIFT[type](n3_2, n3_1, n5_1, n5_2);
        } else if(r > F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = F1_SD_EPS[type](n3_2, n3_1, n5_1, n5_2) * F1_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
        }
    }

    return val;
}

number DNA3Interaction::_f1D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;
    if(r < F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
        if(r > F1_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2 * F1_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * (r - F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
        } else if(r > F1_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            number tmp = exp(-(r - F1_SD_R0[type](n3_2, n3_1, n5_1, n5_2)) * F1_SD_A[type](n3_2, n3_1, n5_1, n5_2));
            val = 2 * (1 - tmp) * tmp * F1_SD_A[type](n3_2, n3_1, n5_1, n5_2);
        } else if(r > F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2 * F1_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * (r - F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
        }
    }

    return F1_SD_EPS[type](n3_2, n3_1, n5_1, n5_2) * val;
}

number DNA3Interaction::_f2_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0.;
    if(r < F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
        if(r > F2_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
            val = F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * F2_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
        } else if(r > F2_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = (F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) / 2.) * (SQR(r - F2_SD_R0[type](n3_2, n3_1, n5_1, n5_2)) - SQR(F2_SD_RC[type](n3_2, n3_1, n5_1, n5_2) - F2_SD_R0[type](n3_2, n3_1, n5_1, n5_2)));
        } else if(r > F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * F2_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
        }
    }
    return val;
}

number DNA3Interaction::_f2D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0.;
    if(r < F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
        if(r > F2_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2. * F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * F2_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * (r - F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
        } else if(r > F2_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * (r - F2_SD_R0[type](n3_2, n3_1, n5_1, n5_2));
        } else if(r > F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2. * F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * F2_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * (r - F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
        }
    }
    return val;
}

number DNA3Interaction::_fakef4_SD(number cost, void *par) {
    if((*(int *)par == CXST_F4_THETA1) && (cost * cost > 1.0001))
        throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", cost);
    if((*(int *)par == CXST_F4_THETA1) && (cost * cost > 1))
        cost = (number)copysign(1, cost);
    number t = Utils::safe_acos(cost);
    return _f4_SD(t, ((int *)par)[0], ((int *)par)[1], ((int *)par)[2], ((int *)par)[3], ((int *)par)[4]);
}

number DNA3Interaction::_fakef4D_SD(number cost, void *par) {
    if((*(int *)par == CXST_F4_THETA1) && (cost * cost > 1.0001))
        throw oxDNAException("In function DNAInteraction::_fakef4() t was found to be out of the range [-1,1] by a large amount, t = %g", cost);
    if((*(int *)par == CXST_F4_THETA1) && (cost * cost > 1))
        cost = (number)copysign(1, cost);
    number t = Utils::safe_acos(cost);
    return -_f4Dsin_SD(t, ((int *)par)[0], ((int *)par)[1], ((int *)par)[2], ((int *)par)[3], ((int *)par)[4]);
}

number DNA3Interaction::_f4_SD(number t, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;
    // std::cout << "f4_SD args: " << t << " " << type << " " << n3 << " " << n5 << std::endl;
    t -= F4_SD_THETA_T0[type](n3_2, n3_1, n5_1, n5_2);
    if(t < 0)
        t *= -1;

    if(t < F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2)) {
        if(t > F4_SD_THETA_TS[type](n3_2, n3_1, n5_1, n5_2)) {
            // smoothing
            val = F4_SD_THETA_B[type](n3_2, n3_1, n5_1, n5_2) * SQR(F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2) - t);
        } 
        else {
            val = (number)1.f - F4_SD_THETA_A[type](n3_2, n3_1, n5_1, n5_2) * SQR(t);
        }
    }

    return val;
}

number DNA3Interaction::_f4D_SD(number t, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;
    number m = (number)1;
    t -= F4_SD_THETA_T0[type](n3_2, n3_1, n5_1, n5_2);
    // this function is a parabola centered in t0. If t < 0 then the value of the function
    // is the same but the value of its derivative has the opposite sign, so m = -1
    if(t < 0) {
        t *= -1;
        m = (number)-1;
    }

    if(t < F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2)) {
        if(t > F4_SD_THETA_TS[type](n3_2, n3_1, n5_1, n5_2)) {
            // smoothing
            val = m * 2 * F4_SD_THETA_B[type](n3_2, n3_1, n5_1, n5_2) * (t - F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2));
        } else
            val = -m * 2 * F4_SD_THETA_A[type](n3_2, n3_1, n5_1, n5_2) * t;
    }

    return val;
}

number DNA3Interaction::_f4Dsin_SD(number t, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;
    number m = (number)1;
    number tt0 = t - F4_SD_THETA_T0[type](n3_2, n3_1, n5_1, n5_2);
    // this function is a parabola centered in t0. If t < 0 then the value of the function
    // is the same but the value of its derivative has the opposite sign, so m = -1
    if(tt0 < 0) {
        tt0 *= -1;
        m = (number)-1;
    }

    if(tt0 < F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2)) {
        number sint = sin(t);
        if(tt0 > F4_SD_THETA_TS[type](n3_2, n3_1, n5_1, n5_2)) {
            // smoothing
            val = m * 2 * F4_SD_THETA_B[type](n3_2, n3_1, n5_1, n5_2) * (tt0 - F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2)) / sint;
        } else {
            if(SQR(sint) > 1e-8)
                val = -m * 2 * F4_SD_THETA_A[type](n3_2, n3_1, n5_1, n5_2) * tt0 / sint;
            else
                val = -m * 2 * F4_SD_THETA_A[type](n3_2, n3_1, n5_1, n5_2);
        }
    }

    return val;
}

number DNA3Interaction::_f5_SD(number f, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;

    if(f > F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2)) {
        if(f < F5_SD_PHI_XS[type](n3_2, n3_1, n5_1, n5_2)) {
            val = F5_SD_PHI_B[type](n3_2, n3_1, n5_1, n5_2) * SQR(F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2) - f);
        } else if(f < 0) {
            val = (number)1.f - F5_SD_PHI_A[type](n3_2, n3_1, n5_1, n5_2) * SQR(f);
        } else
            val = (number)1;
    }

    return val;
}

number DNA3Interaction::_f5D_SD(number f, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    number val = (number)0;

    if(f > F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2)) {
        if(f < F5_SD_PHI_XS[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2 * F5_SD_PHI_B[type](n3_2, n3_1, n5_1, n5_2) * (f - F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2));
        } else if(f < 0) {
            val = (number)-2.f * F5_SD_PHI_A[type](n3_2, n3_1, n5_1, n5_2) * f;
        } else
            val = (number)0;
    }

    return val;
}

void DNA3Interaction::check_input_sanity(std::vector<BaseParticle *> &particles) {
    int N = particles.size();

    for(int i = 0; i < N; i++) {
        BaseParticle *p = particles[i];
        if(p->n3 != P_VIRTUAL && p->n3->index >= N) {
            throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
        }
        if(p->n5 != P_VIRTUAL && p->n5->index >= N) {
            throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
        }

        if(_use_mbf)
            continue;

        // check that the distance between bonded neighbor doesn't exceed a reasonable threshold
        int type_n3_2 = NO_TYPE;
        int type_n5_2 = NO_TYPE;
        if(p->n3 != P_VIRTUAL) {
            BaseParticle *q = p->n3;
            if(q->n3 != P_VIRTUAL) {
                type_n3_2 = q->n3->type;
            }
            if(p->n5 != P_VIRTUAL) {
                type_n5_2 = p->n5->type;
            }

            number mind = _fene_r0_SD(type_n3_2, q->type, p->type, type_n5_2) - _fene_delta_SD(type_n3_2, q->type, p->type, type_n5_2);
            number maxd = _fene_r0_SD(type_n3_2, q->type, p->type, type_n5_2) + _fene_delta_SD(type_n3_2, q->type, p->type, type_n5_2);
            q->set_positions();
            LR_vector rv = p->pos + p->int_centers[DNANucleotide::BACK] - (q->pos + q->int_centers[DNANucleotide::BACK]);
            number r = sqrt(rv * rv);
            if(r > maxd || r < mind) {
                throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf, %lf < d < %lf)", i, p->n3->index, r, mind, maxd);
            }
        }

        type_n3_2 = NO_TYPE;
        type_n5_2 = NO_TYPE;
        if(p->n5 != P_VIRTUAL) {
            BaseParticle *q = p->n5;
            if(p->n3 != P_VIRTUAL) {
                type_n3_2 = p->n3->type;
            }
        if(q->n5 != P_VIRTUAL) {
        type_n5_2 = q->n5->type;
        }

            number mind = _fene_r0_SD(type_n3_2, p->type, q->type, type_n5_2) - _fene_delta_SD(type_n3_2, p->type, q->type, type_n5_2);
            number maxd = _fene_r0_SD(type_n3_2, p->type, q->type, type_n5_2) + _fene_delta_SD(type_n3_2, p->type, q->type, type_n5_2);
            q->set_positions();
            LR_vector rv = p->pos + p->int_centers[DNANucleotide::BACK] - (q->pos + q->int_centers[DNANucleotide::BACK]);
            number r = sqrt(rv * rv);
            if(r > maxd || r < mind) {
                throw oxDNAException("Distance between bonded neighbors %d and %d exceeds acceptable values (d = %lf, %lf < d < %lf)", i, p->n5->index, r, mind, maxd);
            }
        }
    }
}

DNA3Interaction::~DNA3Interaction() {
}

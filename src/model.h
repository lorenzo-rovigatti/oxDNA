/**
 * @file    model.h
 * @date    Jun 2, 2011
 * @author  fromano
 *
 * 
 */

#ifndef MODEL_H_
#define MODEL_H_

/// POSITIONS OF INTERACTION CENTERS
#define POS_BACK -0.4f
#define POS_MM_BACK1 -0.3400f
#define POS_MM_BACK2 0.3408f
#define POS_STACK 0.34f
#define POS_BASE 0.4f
/// POS_STACK - POS_BACK
#define GAMMA 0.74f

/**
 * FENE
 */
#define FENE_EPS 2.0f
#define FENE_R0_OXDNA 0.7525f
#define FENE_R0_OXDNA2 0.7564f
#define FENE_DELTA 0.25f
#define FENE_DELTA2 0.0625f

/**
 * EXCLUDED VOLUME
 */
#define EXCL_EPS 2.0f
#define EXCL_S1 0.70f
#define EXCL_S2 0.33f
#define EXCL_S3 0.515f
#define EXCL_S4 0.515f
#define EXCL_R1 0.675f
#define EXCL_R2 0.32f
#define EXCL_R3 0.50f
#define EXCL_R4 0.50f
#define EXCL_B1 892.016223343f
#define EXCL_B2 4119.70450017f
#define EXCL_B3 1707.30627298f
#define EXCL_B4 1707.30627298f
#define EXCL_RC1 0.711879214356f
#define EXCL_RC2 0.335388426126f
#define EXCL_RC3 0.52329943261f
#define EXCL_RC4 0.52329943261f

/**
 * HYDROGEN BONDING
 */
/// radial part
#define HYDR_F1 0
#define HYDR_EPS_OXDNA 1.077f
#define HYDR_EPS_OXDNA2 1.0678f
#define HYDR_A 8.f
#define HYDR_RC 0.75f
#define HYDR_R0 0.4f
#define HYDR_BLOW -126.243f
#define HYDR_BHIGH -7.87708f
#define HYDR_RLOW 0.34f
#define HYDR_RHIGH 0.7f
#define HYDR_RCLOW 0.276908f
#define HYDR_RCHIGH 0.783775f

/// angular part
#define HYDR_F4_THETA1 2
#define HYDR_F4_THETA2 3
#define HYDR_F4_THETA3 3
#define HYDR_F4_THETA4 4
#define HYDR_F4_THETA7 5
#define HYDR_F4_THETA8 5

#define HYDR_THETA1_A 1.5f
#define HYDR_THETA1_B 4.16038f
#define HYDR_THETA1_T0 0.f
#define HYDR_THETA1_TS 0.7f
#define HYDR_THETA1_TC 0.952381f

#define HYDR_THETA2_A 1.5f
#define HYDR_THETA2_B 4.16038f
#define HYDR_THETA2_T0 0.f
#define HYDR_THETA2_TS 0.7f
#define HYDR_THETA2_TC 0.952381f

#define HYDR_THETA3_A 1.5f
#define HYDR_THETA3_B 4.16038f
#define HYDR_THETA3_T0 0.f
#define HYDR_THETA3_TS 0.7f
#define HYDR_THETA3_TC 0.952381f

#define HYDR_THETA4_A 0.46f
#define HYDR_THETA4_B 0.133855f
#define HYDR_THETA4_T0 PI
#define HYDR_THETA4_TS 0.7f
#define HYDR_THETA4_TC 3.10559f

#define HYDR_THETA7_A 4.f
#define HYDR_THETA7_B 17.0526f
#define HYDR_THETA7_T0 (PI*0.5f)
#define HYDR_THETA7_TS 0.45f
#define HYDR_THETA7_TC 0.555556f

#define HYDR_THETA8_A 4.f
#define HYDR_THETA8_B 17.0526f
#define HYDR_THETA8_T0 (PI*0.5f)
#define HYDR_THETA8_TS 0.45f
#define HYDR_THETA8_TC 0.555556f

/**
 * STACKING
 */
/// radial part
#define STCK_F1 1
#define STCK_BASE_EPS_OXDNA 1.3448f
#define STCK_BASE_EPS_OXDNA2 1.3523f
#define STCK_FACT_EPS_OXDNA 2.6568f
#define STCK_FACT_EPS_OXDNA2 2.6717f
#define STCK_A 6.f
#define STCK_RC 0.9f
#define STCK_R0 0.4f
#define STCK_BLOW -68.1857f
#define STCK_BHIGH -3.12992f
#define STCK_RLOW 0.32f
#define STCK_RHIGH 0.75f
#define STCK_RCLOW 0.23239f
#define STCK_RCHIGH 0.956f

/// angular part, theta
#define STCK_F4_THETA4 0
#define STCK_F4_THETA5 1
#define STCK_F4_THETA6 1
#define STCK_THETA4_A 1.3f
#define STCK_THETA4_B 6.4381f
#define STCK_THETA4_T0 0.f
#define STCK_THETA4_TS 0.8f
#define STCK_THETA4_TC 0.961538f
#define STCK_THETA5_A 0.9f
#define STCK_THETA5_B 3.89361f
#define STCK_THETA5_T0 0.f
#define STCK_THETA5_TS 0.95f
#define STCK_THETA5_TC 1.16959f
#define STCK_THETA6_A 0.9f
#define STCK_THETA6_B 3.89361f
#define STCK_THETA6_T0 0.f
#define STCK_THETA6_TS 0.95f
#define STCK_THETA6_TC 1.16959f

/// angular part, phi
#define STCK_F5_PHI1 0
#define STCK_F5_PHI2 1
#define STCK_PHI1_A 2.0f
#define STCK_PHI1_B 10.9032f
#define STCK_PHI1_XC -0.769231f
#define STCK_PHI1_XS -0.65f
#define STCK_PHI2_A 2.0f
#define STCK_PHI2_B 10.9032f
#define STCK_PHI2_XC -0.769231f
#define STCK_PHI2_XS -0.65f

/**
 *  CROSS STACKING
 */
/// radial part
#define CRST_F2 0
#define CRST_R0 0.575f
#define CRST_RC 0.675f
#define CRST_K 47.5f
#define CRST_BLOW -0.888889f
#define CRST_RLOW 0.495f
#define CRST_RCLOW 0.45f
#define CRST_BHIGH -0.888889f
#define CRST_RHIGH 0.655f
#define CRST_RCHIGH 0.7f
/// angular part;
#define CRST_F4_THETA1 6
#define CRST_F4_THETA2 7
#define CRST_F4_THETA3 7
#define CRST_F4_THETA4 8 
#define CRST_F4_THETA7 9 
#define CRST_F4_THETA8 9 
#define CRST_THETA1_A 2.25f
#define CRST_THETA1_B 7.00545f
#define CRST_THETA1_T0 (PI - 2.35f)
#define CRST_THETA1_TS 0.58f
#define CRST_THETA1_TC 0.766284f
#define CRST_THETA2_A 1.70f
#define CRST_THETA2_B 6.2469f
#define CRST_THETA2_T0 1.f
#define CRST_THETA2_TS 0.68f
#define CRST_THETA2_TC 0.865052f
#define CRST_THETA3_A 1.70f
#define CRST_THETA3_B 6.2469f
#define CRST_THETA3_T0 1.f
#define CRST_THETA3_TS 0.68f
#define CRST_THETA3_TC 0.865052f
#define CRST_THETA4_A 1.50f
#define CRST_THETA4_B 2.59556f
#define CRST_THETA4_T0 0.f
#define CRST_THETA4_TS 0.65f
#define CRST_THETA4_TC 1.02564f
#define CRST_THETA7_A 1.70f
#define CRST_THETA7_B 6.2469f
#define CRST_THETA7_T0 0.875f
#define CRST_THETA7_TS 0.68f
#define CRST_THETA7_TC 0.865052f
#define CRST_THETA8_A 1.70f
#define CRST_THETA8_B 6.2469f
#define CRST_THETA8_T0 0.875f
#define CRST_THETA8_TS 0.68f
#define CRST_THETA8_TC 0.865052f

/**
 *  COAXIAL STACKING
 */
/// radial part
#define CXST_F2 1
#define CXST_R0 0.400f
#define CXST_RC 0.6f
#define CXST_K_OXDNA 46.0f
#define CXST_K_OXDNA2 58.5f
#define CXST_BLOW -2.13158f
#define CXST_RLOW 0.22f
#define CXST_RCLOW 0.177778f
#define CXST_BHIGH -2.13158f
#define CXST_RHIGH 0.58f
#define CXST_RCHIGH 0.6222222f
/// angular part;
#define CXST_F4_THETA1 10
#define CXST_F4_THETA4 11
#define CXST_F4_THETA5 12 
#define CXST_F4_THETA6 12 
#define CXST_THETA1_A   2.f
#define CXST_THETA1_B   10.9032f
#define CXST_THETA1_T0_OXDNA (PI - 0.60f)
#define CXST_THETA1_T0_OXDNA2 (PI - 0.25f)
#define CXST_THETA1_TS  0.65f
#define CXST_THETA1_TC  0.769231f
#define CXST_THETA1_SA  20.f
#define CXST_THETA1_SB (PI - 0.1f*(PI - (PI - 0.25f)))
#define CXST_THETA4_A   1.3f
#define CXST_THETA4_B   6.4381f
#define CXST_THETA4_T0  0.f
#define CXST_THETA4_TS  0.8f
#define CXST_THETA4_TC  0.961538f
#define CXST_THETA5_A   0.9f
#define CXST_THETA5_B 	3.89361f  
#define CXST_THETA5_T0  0.f
#define CXST_THETA5_TS  0.95f
#define CXST_THETA5_TC  1.16959f
#define CXST_THETA6_A   0.9f
#define CXST_THETA6_B 	3.89361f  
#define CXST_THETA6_T0  0.f
#define CXST_THETA6_TS  0.95f
#define CXST_THETA6_TC  1.16959f
/// angular part, phi3 and phi4
#define CXST_F5_PHI3 2
#define CXST_F5_PHI4 3
#define CXST_PHI3_A 2.0f
#define CXST_PHI3_B 10.9032f
#define CXST_PHI3_XC -0.769231f
#define CXST_PHI3_XS -0.65f
#define CXST_PHI4_A 2.0f
#define CXST_PHI4_B 10.9032f
#define CXST_PHI4_XC -0.769231f
#define CXST_PHI4_XS -0.65f

#define HYDR_T1_MESH_POINTS 6   // perfect
#define HYDR_T2_MESH_POINTS 6   // perfect
#define HYDR_T3_MESH_POINTS HYDR_T2_MESH_POINTS
#define HYDR_T4_MESH_POINTS 50  // almost perfect 
#define HYDR_T7_MESH_POINTS 12  // almost perfect
#define HYDR_T8_MESH_POINTS HYDR_T7_MESH_POINTS

#define STCK_T4_MESH_POINTS 6   // perfect
#define STCK_T5_MESH_POINTS 6   // perfect 
#define STCK_T6_MESH_POINTS STCK_T5_MESH_POINTS

#define CRST_T1_MESH_POINTS 250 // good enough
#define CRST_T2_MESH_POINTS 80  // almost perfect
#define CRST_T3_MESH_POINTS CRST_T2_MESH_POINTS
#define CRST_T4_MESH_POINTS 6   // perfect
#define CRST_T7_MESH_POINTS 250 // good enough
#define CRST_T8_MESH_POINTS CRST_T7_MESH_POINTS

#define CXST_T1_MESH_POINTS 250 // perfetto
#define CXST_T4_MESH_POINTS 6   // perfetto
#define CXST_T5_MESH_POINTS 6   // perfetto
#define CXST_T6_MESH_POINTS CXST_T5_MESH_POINTS


#endif /* MODEL_H_ */
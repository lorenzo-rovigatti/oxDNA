/**
 * 
 */

#ifndef HYBRID_MODEL_H_
#define HYBRID_MODEL_H_



//------------------PARAMETERS FOR HYBRID INTERACTION------------------//
/**
 * HYDROGEN BONDING
 */
/// radial part
#define HYDR_F1_hybrid 0
#define HYDR_EPS_hybrid 1.43f
#define HYDR_A_hybrid 8.f
#define HYDR_RC_hybrid 0.75f
#define HYDR_R0_hybrid 0.4f
#define HYDR_BLOW_hybrid -126.243f
#define HYDR_BHIGH_hybrid -7.87708f
#define HYDR_RLOW_hybrid 0.34f
#define HYDR_RHIGH_hybrid 0.7f
#define HYDR_RCLOW_hybrid 0.276908f
#define HYDR_RCHIGH_hybrid 0.783775f

/// angular part
#define HYDR_F4_THETA1_hybrid 2
#define HYDR_F4_THETA2_hybrid 3
#define HYDR_F4_THETA3_hybrid 3
#define HYDR_F4_THETA4_hybrid 4
#define HYDR_F4_THETA7_hybrid 5
#define HYDR_F4_THETA8_hybrid 5

#define HYDR_THETA1_A_hybrid 1.5f
#define HYDR_THETA1_B_hybrid 4.16038f
#define HYDR_THETA1_T0_hybrid 0.f
#define HYDR_THETA1_TS_hybrid 0.7f
#define HYDR_THETA1_TC_hybrid 0.952381f

#define HYDR_THETA2_A_hybrid 1.5f
#define HYDR_THETA2_B_hybrid 4.16038f
#define HYDR_THETA2_T0_hybrid 0.f
#define HYDR_THETA2_TS_hybrid 0.7f
#define HYDR_THETA2_TC_hybrid 0.952381f

#define HYDR_THETA3_A_hybrid 1.5f
#define HYDR_THETA3_B_hybrid 4.16038f
#define HYDR_THETA3_T0_hybrid 0.f
#define HYDR_THETA3_TS_hybrid 0.7f
#define HYDR_THETA3_TC_hybrid 0.952381f

#define HYDR_THETA4_A_hybrid 0.46f
#define HYDR_THETA4_B_hybrid 0.133855f
#define HYDR_THETA4_T0_hybrid PI
#define HYDR_THETA4_TS_hybrid 0.7f
#define HYDR_THETA4_TC_hybrid 3.10559f

#define HYDR_THETA7_A_hybrid 4.f
#define HYDR_THETA7_B_hybrid 17.0526f
#define HYDR_THETA7_T0_hybrid (PI*0.5f)
#define HYDR_THETA7_TS_hybrid 0.45f
#define HYDR_THETA7_TC_hybrid 0.555556f

#define HYDR_THETA8_A_hybrid 4.f
#define HYDR_THETA8_B_hybrid 17.0526f
#define HYDR_THETA8_T0_hybrid (PI*0.5f)
#define HYDR_THETA8_TS_hybrid 0.45f
#define HYDR_THETA8_TC_hybrid 0.555556f


/**
 *  CROSS STACKING
 */
/// radial part
#define CRST_F2_hybrid 0  
#define CRST_R0_hybrid 0.575f
#define CRST_RC_hybrid 0.675f
#define CRST_K_hybrid 44.535f
#define CRST_BLOW_hybrid -0.888889f
#define CRST_RLOW_hybrid 0.495f
#define CRST_RCLOW_hybrid 0.45f
#define CRST_BHIGH_hybrid -0.888889f
#define CRST_RHIGH_hybrid 0.655f
#define CRST_RCHIGH_hybrid 0.7f
/// angular part;
#define CRST_F4_THETA1_hybrid 6
#define CRST_F4_THETA2_hybrid 7
#define CRST_F4_THETA3_hybrid 7
#define CRST_F4_THETA4_hybrid 8 
#define CRST_F4_THETA7_hybrid 9 
#define CRST_F4_THETA8_hybrid 9  
#define CRST_THETA1_A_hybrid 2.25f
#define CRST_THETA1_B_hybrid 7.00545f
#define CRST_THETA1_T0_hybrid (PI - 2.35f)
#define CRST_THETA1_TS_hybrid 0.58f
#define CRST_THETA1_TC_hybrid 0.766284f
#define CRST_THETA2_A_hybrid 1.70f
#define CRST_THETA2_B_hybrid 6.2469f
#define CRST_THETA2_T0_hybrid 1.f
#define CRST_THETA2_TS_hybrid 0.68f
#define CRST_THETA2_TC_hybrid 0.865052f
#define CRST_THETA3_A_hybrid 1.70f
#define CRST_THETA3_B_hybrid 6.2469f
#define CRST_THETA3_T0_hybrid 1.f
#define CRST_THETA3_TS_hybrid 0.68f
#define CRST_THETA3_TC_hybrid 0.865052f
#define CRST_THETA4_A_hybrid 1.50f
#define CRST_THETA4_B_hybrid 2.59556f
#define CRST_THETA4_T0_hybrid 0.f
#define CRST_THETA4_TS_hybrid 0.65f
#define CRST_THETA4_TC_hybrid 1.02564f
#define CRST_THETA7_A_hybrid 1.70f
#define CRST_THETA7_B_hybrid 6.2469f
#define CRST_THETA7_T0_hybrid 0.875f
#define CRST_THETA7_TS_hybrid 0.68f
#define CRST_THETA7_TC_hybrid 0.865052f
#define CRST_THETA8_A_hybrid 1.70f
#define CRST_THETA8_B_hybrid 6.2469f
#define CRST_THETA8_T0_hybrid 0.875f
#define CRST_THETA8_TS_hybrid 0.68f
#define CRST_THETA8_TC_hybrid 0.865052f

//Mesh stuff
#define HYDR_T1_MESH_POINTS_hybrid 6   // perfect
#define HYDR_T2_MESH_POINTS_hybrid 6   // perfect
#define HYDR_T3_MESH_POINTS_hybrid HYDR_T2_MESH_POINTS_hybrid
#define HYDR_T4_MESH_POINTS_hybrid 50  // almost perfect 
#define HYDR_T7_MESH_POINTS_hybrid 12  // almost perfect
#define HYDR_T8_MESH_POINTS_hybrid HYDR_T7_MESH_POINTS_hybrid

#define CRST_T1_MESH_POINTS_hybrid 250 // good enough
#define CRST_T2_MESH_POINTS_hybrid 80  // almost perfect
#define CRST_T3_MESH_POINTS_hybrid CRST_T2_MESH_POINTS_hybrid
#define CRST_T4_MESH_POINTS_hybrid 6   // perfect
#define CRST_T7_MESH_POINTS_hybrid 250 // good enough
#define CRST_T8_MESH_POINTS_hybrid CRST_T7_MESH_POINTS_hybrid

#endif 

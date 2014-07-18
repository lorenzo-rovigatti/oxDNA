/*
 * interaction_consts.h
 *
 */

#ifndef RNAMODEL_H_
#define RNAMODEL_H_

#include <cstdlib>
#include <cstdio>

#include "../Utilities/parse_input/parse_input.h"

//#define PI 3.141592653589793238462643f

#define RNA_MODEL_VERSION "DEC2013"

#define RNA_HYDR_F1 0
#define RNA_HYDR_F4_THETA1 2
#define RNA_HYDR_F4_THETA2 3
#define RNA_HYDR_F4_THETA3 3
#define RNA_HYDR_F4_THETA4 4
#define RNA_HYDR_F4_THETA7 5
#define RNA_HYDR_F4_THETA8 5
#define RNA_STCK_F1 1

#define RNA_STCK_F4_THETA4 0
#define RNA_STCK_F4_THETA5 1
#define RNA_STCK_F4_THETA6 13

#define RNA_STCK_F4_THETAB1 14
#define RNA_STCK_F4_THETAB2 15

#define RNA_STCK_F5_PHI1 0
#define RNA_STCK_F5_PHI2 1

#define RNA_CRST_F2 0

#define RNA_CRST_F4_THETA1 6
#define RNA_CRST_F4_THETA2 7
#define RNA_CRST_F4_THETA3 7
#define RNA_CRST_F4_THETA4 8
#define RNA_CRST_F4_THETA7 9
#define RNA_CRST_F4_THETA8 9

#define RNA_CXST_F2 1

#define RNA_CXST_F4_THETA1 10
#define RNA_CXST_F4_THETA4 11
#define RNA_CXST_F4_THETA5 12
#define RNA_CXST_F4_THETA6 12

#define RNA_CXST_F5_PHI3 2
#define RNA_CXST_F5_PHI4 3

struct Model {

	float RNA_POS_BACK;
	float RNA_POS_STACK;
	float RNA_POS_BASE;
	// RNA_POS_STACK - RNA_POS_BACK
	float RNA_GAMMA;

	//this is stacking site with your 3-neighbor
	float	RNA_POS_STACK_3_a1;
	float	RNA_POS_STACK_3_a2;
	//this is stacking site for interaction with your 5-neighbor
	float	RNA_POS_STACK_5_a1 ;
	float	RNA_POS_STACK_5_a2 ;



	/*
	 * FENE
	 */
	float RNA_FENE_EPS;
	float RNA_FENE_R0;
	float RNA_FENE_DELTA;
	float RNA_FENE_DELTA2;

	/*
	 * EXCLUDED VOLUME
	 */
	float RNA_EXCL_EPS;
	float RNA_EXCL_S1;
	float RNA_EXCL_S2;
	float RNA_EXCL_S3;
	float RNA_EXCL_S4;
	float RNA_EXCL_R1;
	float RNA_EXCL_R2;
	float RNA_EXCL_R3;
	float RNA_EXCL_R4;
	float RNA_EXCL_B1;
	float RNA_EXCL_B2;
	float RNA_EXCL_B3;
	float RNA_EXCL_B4;
	float RNA_EXCL_RC1;
	float RNA_EXCL_RC2;
	float RNA_EXCL_RC3;
	float RNA_EXCL_RC4;

	/*
	 * HYDROGEN BONDING
	 */
	// radial part
	float RNA_HYDR_EPS;
	float RNA_HYDR_A;
	float RNA_HYDR_RC;
	float RNA_HYDR_R0;
	float RNA_HYDR_BLOW;
	float RNA_HYDR_BHIGH;
	float RNA_HYDR_RLOW;
	float RNA_HYDR_RHIGH;
	float RNA_HYDR_RCLOW;
	float RNA_HYDR_RCHIGH;

	// angular part

	float RNA_HYDR_THETA1_A;
	float RNA_HYDR_THETA1_B;
	float RNA_HYDR_THETA1_T0;
	float RNA_HYDR_THETA1_TS;
	float RNA_HYDR_THETA1_TC;

	float RNA_HYDR_THETA2_A;
	float RNA_HYDR_THETA2_B;
	float RNA_HYDR_THETA2_T0;
	float RNA_HYDR_THETA2_TS;
	float RNA_HYDR_THETA2_TC;

	float RNA_HYDR_THETA3_A;
	float RNA_HYDR_THETA3_B;
	float RNA_HYDR_THETA3_T0;
	float RNA_HYDR_THETA3_TS;
	float RNA_HYDR_THETA3_TC;

	float RNA_HYDR_THETA4_A;
	float RNA_HYDR_THETA4_B;
	float RNA_HYDR_THETA4_T0;
	float RNA_HYDR_THETA4_TS;
	float RNA_HYDR_THETA4_TC;

	float RNA_HYDR_THETA7_A;
	float RNA_HYDR_THETA7_B;
	float RNA_HYDR_THETA7_T0;
	float RNA_HYDR_THETA7_TS;
	float RNA_HYDR_THETA7_TC;

	float RNA_HYDR_THETA8_A;
	float RNA_HYDR_THETA8_B;
	float RNA_HYDR_THETA8_T0;
	float RNA_HYDR_THETA8_TS;
	float RNA_HYDR_THETA8_TC;

	/*
	 * STACKING
	 */
	// radial part
	float RNA_STCK_BASE_EPS;
	float RNA_STCK_FACT_EPS;
	float RNA_STCK_A;
	float RNA_STCK_RC;
	float RNA_STCK_R0;
	float RNA_STCK_BLOW;
	float RNA_STCK_BHIGH;
	float RNA_STCK_RLOW;
	float RNA_STCK_RHIGH;
	float RNA_STCK_RCLOW;
	float RNA_STCK_RCHIGH;

	// angular part, theta
	float RNA_STCK_THETA4_A;
	float RNA_STCK_THETA4_B;
	float RNA_STCK_THETA4_T0;
	float RNA_STCK_THETA4_TS;
	float RNA_STCK_THETA4_TC;
	float RNA_STCK_THETA5_A;
	float RNA_STCK_THETA5_B;
	float RNA_STCK_THETA5_T0;
	float RNA_STCK_THETA5_TS;
	float RNA_STCK_THETA5_TC;
	float RNA_STCK_THETA6_A;
	float RNA_STCK_THETA6_B;
	float RNA_STCK_THETA6_T0;
	float RNA_STCK_THETA6_TS;
	float RNA_STCK_THETA6_TC;

	float STCK_THETAB1_A;
	float STCK_THETAB1_B;
	float STCK_THETAB1_T0;
	float STCK_THETAB1_TS;
	float STCK_THETAB1_TC;

	float STCK_THETAB2_A;
	float STCK_THETAB2_B;
	float STCK_THETAB2_T0;
	float STCK_THETAB2_TS;
	float STCK_THETAB2_TC;




	// angular part, phi
	float RNA_STCK_PHI1_A;
	float RNA_STCK_PHI1_B;
	float RNA_STCK_PHI1_XC;
	float RNA_STCK_PHI1_XS;
	float RNA_STCK_PHI2_A;
	float RNA_STCK_PHI2_B;
	float RNA_STCK_PHI2_XC;
	float RNA_STCK_PHI2_XS;

	/*
	 *  CROSS STACKING
	 */
	// radial part
	float RNA_CRST_R0;
	float RNA_CRST_RC;
	float RNA_CRST_K;
	float RNA_CRST_BLOW;
	float RNA_CRST_RLOW;
	float RNA_CRST_RCLOW;
	float RNA_CRST_BHIGH;
	float RNA_CRST_RHIGH;
	float RNA_CRST_RCHIGH;
	// angular part;
	float RNA_CRST_THETA1_A;
	float RNA_CRST_THETA1_B;
	float RNA_CRST_THETA1_T0;
	float RNA_CRST_THETA1_TS;
	float RNA_CRST_THETA1_TC;
	float RNA_CRST_THETA2_A;
	float RNA_CRST_THETA2_B;
	float RNA_CRST_THETA2_T0;
	float RNA_CRST_THETA2_TS;
	float RNA_CRST_THETA2_TC;
	float RNA_CRST_THETA3_A;
	float RNA_CRST_THETA3_B;
	float RNA_CRST_THETA3_T0;
	float RNA_CRST_THETA3_TS;
	float RNA_CRST_THETA3_TC;
	float RNA_CRST_THETA4_A;
	float RNA_CRST_THETA4_B;
	float RNA_CRST_THETA4_T0;
	float RNA_CRST_THETA4_TS;
	float RNA_CRST_THETA4_TC;
	float RNA_CRST_THETA7_A;
	float RNA_CRST_THETA7_B;
	float RNA_CRST_THETA7_T0;
	float RNA_CRST_THETA7_TS;
	float RNA_CRST_THETA7_TC;
	float RNA_CRST_THETA8_A;
	float RNA_CRST_THETA8_B;
	float RNA_CRST_THETA8_T0;
	float RNA_CRST_THETA8_TS;
	float RNA_CRST_THETA8_TC;

	/*
	 *  COAXIAL STACKING
	 */
	// radial part
	float RNA_CXST_R0;
	float RNA_CXST_RC;
	float RNA_CXST_K;
	float RNA_CXST_BLOW;
	float RNA_CXST_RLOW;
	float RNA_CXST_RCLOW;
	float RNA_CXST_BHIGH;
	float RNA_CXST_RHIGH;
	float RNA_CXST_RCHIGH;
	// angular part;

	float RNA_CXST_THETA1_A;
	float RNA_CXST_THETA1_B;
	float RNA_CXST_THETA1_T0;
	float RNA_CXST_THETA1_TS;
	float RNA_CXST_THETA1_TC;
	float RNA_CXST_THETA4_A;
	float RNA_CXST_THETA4_B;
	float RNA_CXST_THETA4_T0;
	float RNA_CXST_THETA4_TS;
	float RNA_CXST_THETA4_TC;
	float RNA_CXST_THETA5_A;
	float RNA_CXST_THETA5_B;
	float RNA_CXST_THETA5_T0;
	float RNA_CXST_THETA5_TS;
	float RNA_CXST_THETA5_TC;
	float RNA_CXST_THETA6_A;
	float RNA_CXST_THETA6_B;
	float RNA_CXST_THETA6_T0;
	float RNA_CXST_THETA6_TS;
	float RNA_CXST_THETA6_TC;
	// angular part, phi3 and phi4

	float RNA_CXST_PHI3_A;
	float RNA_CXST_PHI3_B;
	float RNA_CXST_PHI3_XC;
	float RNA_CXST_PHI3_XS;
	float RNA_CXST_PHI4_A;
	float RNA_CXST_PHI4_B;
	float RNA_CXST_PHI4_XC;
	float RNA_CXST_PHI4_XS;

//new stuff
	float INCLINATION;
    float BASE_SIZE_R;
    float ROT_SUGAR;
    float p3_x, p3_y, p3_z, p5_x, p5_y, p5_z;
    float RNA_POS_BACK_a1, RNA_POS_BACK_a2, RNA_POS_BACK_a3;

	Model() {
// POSITIONS OF INTERACTION CENTERS
		RNA_POS_BACK = -0.4f;
		RNA_POS_STACK = 0.34f;
		RNA_POS_BASE = 0.4f;


//NEW STUFF
		//this is stacking site with you 3-neighbor
		//RNA_POS_STACK_3_a1 = 0.34f;
	   //	RNA_POS_STACK_3_a2 = 0.0f;

		//this is stacking site for interaction with you 5-neighbor
		//RNA_POS_STACK_5_a1 = 0.34f;
		//RNA_POS_STACK_5_a2 = 0.0f;


		RNA_POS_STACK_3_a1 =  0.4 ;
		RNA_POS_STACK_3_a2 =  0.1 ;
		RNA_POS_STACK_5_a1 =   0.124906078525 ;
		RNA_POS_STACK_5_a2 =  -0.00866274917473;

		RNA_POS_BACK_a1 = -0.4;
		RNA_POS_BACK_a3 = 0.2;
		RNA_POS_BACK_a2 = 0.0;

		//RNA_POS_BACK_a1 = -0.4f;
		//RNA_POS_BACK_a2 = 0;
		//RNA_POS_BACK_a3 = 0;

		INCLINATION = 0.0f;
		BASE_SIZE_R = 0.01f;
		ROT_SUGAR = 0.f; //1.169f;
		//p3_x=0, p3_y=0, p3_z=1, p5_x=0, p5_y=0, p5_z=1;

		p5_x = -0.104402;
		p5_y = -0.841783 ;
		p5_z = 0.529624 ;
		p3_x = -0.462510 ;
		p3_y = -0.528218 ;
		p3_z = 0.712089 ;

// RNA_POS_STACK - RNA_POS_BACK
		RNA_GAMMA = 0.74f;

		/*
		 * FENE
		 */
		RNA_FENE_EPS = 2.0f;
		//RNA_FENE_R0 = 0.7525f;
		RNA_FENE_R0 =  0.761070781051 ;
		RNA_FENE_DELTA = 0.25f;
		RNA_FENE_DELTA2 = 0.0625f;

		/*
		 * EXCLUDED VOLUME
		 */
		RNA_EXCL_EPS = 2.0f;
		RNA_EXCL_S1 = 0.70f;
		RNA_EXCL_S2 = 0.33f;
		RNA_EXCL_S3 = 0.515f;
		RNA_EXCL_S4 = 0.515f;
		RNA_EXCL_R1 = 0.675f;
		RNA_EXCL_R2 = 0.32f;
		RNA_EXCL_R3 = 0.50f;
		RNA_EXCL_R4 = 0.50f;
		RNA_EXCL_B1 = 892.016223343f;
		RNA_EXCL_B2 = 4119.70450017f;
		RNA_EXCL_B3 = 1707.30627298f;
		RNA_EXCL_B4 = 1707.30627298f;
		RNA_EXCL_RC1 = 0.711879214356f;
		RNA_EXCL_RC2 = 0.335388426126f;
		RNA_EXCL_RC3 = 0.52329943261f;
		RNA_EXCL_RC4 = 0.52329943261f;

		/*
		 * HYDROGEN BONDING
		 */
// radial part
		//RNA_HYDR_EPS = 1.077f;
		RNA_HYDR_EPS = 0.870439f;
		RNA_HYDR_A = 8.f;
		RNA_HYDR_RC = 0.75f;
		RNA_HYDR_R0 = 0.4f;
		RNA_HYDR_BLOW = -126.243f;
		RNA_HYDR_BHIGH = -7.87708f;
		RNA_HYDR_RLOW = 0.34f;
		RNA_HYDR_RHIGH = 0.7f;
		RNA_HYDR_RCLOW = 0.276908f;
		RNA_HYDR_RCHIGH = 0.783775f;

// angular part

		RNA_HYDR_THETA1_A = 1.5f;
		RNA_HYDR_THETA1_B = 4.16038f;
		RNA_HYDR_THETA1_T0 = 0.f;
		RNA_HYDR_THETA1_TS = 0.7f;
		RNA_HYDR_THETA1_TC = 0.952381f;

		RNA_HYDR_THETA2_A = 1.5f;
		RNA_HYDR_THETA2_B = 4.16038f;
		RNA_HYDR_THETA2_T0 = 0.f;
		RNA_HYDR_THETA2_TS = 0.7f;
		RNA_HYDR_THETA2_TC = 0.952381f;

		RNA_HYDR_THETA3_A = 1.5f;
		RNA_HYDR_THETA3_B = 4.16038f;
		RNA_HYDR_THETA3_T0 = 0.f;
		RNA_HYDR_THETA3_TS = 0.7f;
		RNA_HYDR_THETA3_TC = 0.952381f;

		RNA_HYDR_THETA4_A = 0.46f;
		RNA_HYDR_THETA4_B = 0.133855f;
		RNA_HYDR_THETA4_T0 = PI;
		RNA_HYDR_THETA4_TS = 0.7f;
		RNA_HYDR_THETA4_TC = 3.10559f;

		RNA_HYDR_THETA7_A = 4.f;
		RNA_HYDR_THETA7_B = 17.0526f;
		RNA_HYDR_THETA7_T0 = (PI * 0.5f);
		RNA_HYDR_THETA7_TS = 0.45f;
		RNA_HYDR_THETA7_TC = 0.555556f;

		RNA_HYDR_THETA8_A = 4.f;
		RNA_HYDR_THETA8_B = 17.0526f;
		RNA_HYDR_THETA8_T0 = (PI * 0.5f);
		RNA_HYDR_THETA8_TS = 0.45f;
		RNA_HYDR_THETA8_TC = 0.555556f;

		/*
		 * STACKING
		 */
// radial part
		//RNA_STCK_BASE_EPS = 1.3448f;
		//RNA_STCK_FACT_EPS = 2.6568f;
		RNA_STCK_BASE_EPS = 1.40206f;
		RNA_STCK_FACT_EPS = 2.77f;

		RNA_STCK_RC = 0.93;
		RNA_STCK_RLOW = 0.35;
		RNA_STCK_RCLOW = 0.26239;
		RNA_STCK_RHIGH = 0.78;
		RNA_STCK_RCHIGH = 0.986;
		RNA_STCK_R0 =  0.43;

		RNA_STCK_A = 6.f;
		//RNA_STCK_RC = 0.9f;
		//RNA_STCK_R0 = 0.4f;
		RNA_STCK_BLOW = -68.1857f;
		RNA_STCK_BHIGH = -3.12992f;
	//	RNA_STCK_RLOW = 0.32f;
	//	RNA_STCK_RHIGH = 0.75f;
	//	RNA_STCK_RCLOW = 0.23239f;
	//	RNA_STCK_RCHIGH = 0.956f;

// angular part, theta
		RNA_STCK_THETA4_A = 0;
		RNA_STCK_THETA4_B = 0;

		//RNA_STCK_THETA4_A = 1.3f;
		//RNA_STCK_THETA4_B = 6.4381f;
		RNA_STCK_THETA4_T0 = 0.f;
		RNA_STCK_THETA4_TS = 0.8f;
		RNA_STCK_THETA4_TC = 0.961538f;
		RNA_STCK_THETA5_A = 0.9f;
		RNA_STCK_THETA5_B = 3.89361f;
		RNA_STCK_THETA5_T0 = 0.f;
		RNA_STCK_THETA5_TS = 0.95f;
		RNA_STCK_THETA5_TC = 1.16959f;
		RNA_STCK_THETA6_A = 0.9f;
		RNA_STCK_THETA6_B = 3.89361f;
		RNA_STCK_THETA6_T0 = 0.f;
		RNA_STCK_THETA6_TS = 0.95f;
		RNA_STCK_THETA6_TC = 1.16959f;


		STCK_THETAB2_A = 1.3;
		STCK_THETAB2_B = 6.4381;
		STCK_THETAB2_TS = 0.8;
		STCK_THETAB2_TC = 0.961538;
		STCK_THETAB1_A = 1.3;
		STCK_THETAB1_B = 6.4381;
		STCK_THETAB1_TS = 0.8;
		STCK_THETAB1_TC = 0.961538;

		//STCK_THETAB1_A = 0.9f;
		//STCK_THETAB1_B = 3.89361f;
		STCK_THETAB1_T0 = 0.f;
		//STCK_THETAB1_TS = 0.95f;
		//STCK_THETAB1_TC = 1.16959f;
		//STCK_THETAB1_A = 0.9f;
		//STCK_THETAB1_B = 3.89361f;
		//STCK_THETAB1_T0 = 0.f;
		//STCK_THETAB1_TS = 0.95f;
		//STCK_THETAB1_TC = 1.16959f;

		//STCK_THETAB1_A = 0.9f;
		//STCK_THETAB1_B = 3.89361f;
		//STCK_THETAB1_T0 = 0.f;
		//STCK_THETAB1_TS = 0.95f;
		//STCK_THETAB1_TC = 1.16959f;
		//STCK_THETAB2_A = 0.9f;
		//STCK_THETAB2_B = 3.89361f;
		STCK_THETAB2_T0 = 0.f;
		//STCK_THETAB2_TS = 0.95f;
		//STCK_THETAB2_TC = 1.16959f;



// angular part, phi
		RNA_STCK_PHI1_A = 2.0f;
		RNA_STCK_PHI1_B = 10.9032f;
		RNA_STCK_PHI1_XC = -0.769231f;
		RNA_STCK_PHI1_XS = -0.65f;
		RNA_STCK_PHI2_A = 2.0f;
		RNA_STCK_PHI2_B = 10.9032f;
		RNA_STCK_PHI2_XC = -0.769231f;
		RNA_STCK_PHI2_XS = -0.65f;

		/*
		 *  CROSS STACKING
		 */
// radial part
		//RNA_CRST_R0 = 0.575f;
		//RNA_CRST_RC = 0.675f;
		//RNA_CRST_K = 47.5f;
		RNA_CRST_BLOW = -0.888889f;
		//RNA_CRST_RLOW = 0.495f;
		//RNA_CRST_RCLOW = 0.45f;
		RNA_CRST_BHIGH = -0.888889f;
		//RNA_CRST_RHIGH = 0.655f;
		//RNA_CRST_RCHIGH = 0.7f;

		RNA_CRST_RC = 0.6;
		RNA_CRST_RLOW = 0.42;
		RNA_CRST_RCLOW = 0.375;
		RNA_CRST_RHIGH = 0.58;
		RNA_CRST_RCHIGH = 0.625;
		RNA_CRST_R0 =  0.5;
		RNA_CRST_K = 59.9626;
		RNA_CRST_THETA1_T0 = 0.505;
		RNA_CRST_THETA2_T0 = 1.266;
		RNA_CRST_THETA3_T0 = 1.266;
		RNA_CRST_THETA4_T0 = 0;
		RNA_CRST_THETA7_T0 = 0.309;
		RNA_CRST_THETA8_T0 = 0.309;
// angular part;
		RNA_CRST_THETA1_A = 2.25f;
		RNA_CRST_THETA1_B = 7.00545f;
		//RNA_CRST_THETA1_T0 = (PI - 2.35f);
		RNA_CRST_THETA1_TS = 0.58f;
		RNA_CRST_THETA1_TC = 0.766284f;
		RNA_CRST_THETA2_A = 1.70f;
		RNA_CRST_THETA2_B = 6.2469f;
		//RNA_CRST_THETA2_T0 = 1.f;
		RNA_CRST_THETA2_TS = 0.68f;
		RNA_CRST_THETA2_TC = 0.865052f;
		RNA_CRST_THETA3_A = 1.70f;
		RNA_CRST_THETA3_B = 6.2469f;
		//RNA_CRST_THETA3_T0 = 1.f;
		RNA_CRST_THETA3_TS = 0.68f;
		RNA_CRST_THETA3_TC = 0.865052f;
		RNA_CRST_THETA4_A = 1.50f;
		RNA_CRST_THETA4_B = 2.59556f;
		//RNA_CRST_THETA4_T0 = 0.f;
		RNA_CRST_THETA4_TS = 0.65f;
		RNA_CRST_THETA4_TC = 1.02564f;
		RNA_CRST_THETA7_A = 1.70f;
		RNA_CRST_THETA7_B = 6.2469f;
		//RNA_CRST_THETA7_T0 = 0.875f;
		RNA_CRST_THETA7_TS = 0.68f;
		RNA_CRST_THETA7_TC = 0.865052f;
		RNA_CRST_THETA8_A = 1.70f;
		RNA_CRST_THETA8_B = 6.2469f;
		//RNA_CRST_THETA8_T0 = 0.875f;
		RNA_CRST_THETA8_TS = 0.68f;
		RNA_CRST_THETA8_TC = 0.865052f;

		/*
		 *  COAXIAL STACKING
		 */
// radial part


		RNA_CXST_K = 80.;
		RNA_CXST_RC = 0.6;
		RNA_CXST_RLOW = 0.42;
		RNA_CXST_RCLOW = 0.375;
		RNA_CXST_RHIGH = 0.58;
		RNA_CXST_RCHIGH = 0.625;
		RNA_CXST_R0 =  0.5;
		RNA_CXST_BLOW = -0.888889;
		RNA_CXST_BHIGH = -0.888889;



		//RNA_CXST_R0 = 0.400f;
		//RNA_CXST_RC = 0.6f;
		//RNA_CXST_K = 46.0f;
		//RNA_CXST_BLOW = -2.13158f;
		//RNA_CXST_RLOW = 0.22f;
		//RNA_CXST_RCLOW = 0.177778f;
		//RNA_CXST_BHIGH = -2.13158f;
		//RNA_CXST_RHIGH = 0.58f;
		//RNA_CXST_RCHIGH = 0.6222222f;
// angular part;




		RNA_CXST_THETA1_T0 = 2.592;
		RNA_CXST_THETA4_T0 = 0.151;
		RNA_CXST_THETA5_T0 = 0.685;
		RNA_CXST_THETA6_T0 = 0.685;//2.523;


		RNA_CXST_THETA1_A = 2.f;
		RNA_CXST_THETA1_B = 10.9032f;
		//RNA_CXST_THETA1_T0 = (PI - 0.60f);
		RNA_CXST_THETA1_TS = 0.65f;
		RNA_CXST_THETA1_TC = 0.769231f;
		RNA_CXST_THETA4_A = 1.3f;
		RNA_CXST_THETA4_B = 6.4381f;
		//RNA_CXST_THETA4_T0 = 0.f;
		RNA_CXST_THETA4_TS = 0.8f;
		RNA_CXST_THETA4_TC = 0.961538f;
		RNA_CXST_THETA5_A = 0.9f;
		RNA_CXST_THETA5_B = 3.89361f;
		//RNA_CXST_THETA5_T0 = 0.f;
		RNA_CXST_THETA5_TS = 0.95f;
		RNA_CXST_THETA5_TC = 1.16959f;
		RNA_CXST_THETA6_A = 0.9f;
		RNA_CXST_THETA6_B = 3.89361f;
		//RNA_CXST_THETA6_T0 = 0.f;
		RNA_CXST_THETA6_TS = 0.95f;
		RNA_CXST_THETA6_TC = 1.16959f;
// angular part, phi3 and phi4

		RNA_CXST_PHI3_A = 2.0f;
		RNA_CXST_PHI3_B = 10.9032f;
		RNA_CXST_PHI3_XC = -0.769231f;
		RNA_CXST_PHI3_XS = -0.65f;
		RNA_CXST_PHI4_A = 2.0f;
		RNA_CXST_PHI4_B = 10.9032f;
		RNA_CXST_PHI4_XC = -0.769231f;
		RNA_CXST_PHI4_XS = -0.65f;

	}
	void load_model_from_file(const char *fname) {
		FILE *f = fopen(fname, "r");
		if (f == NULL) {
			fprintf(stderr, "Cannot open %s", fname);
			return;
		}

		input_file input;
		//char type_str[512];
		//char name_str[512];
		float tmp_value;
	//	int type;
		loadInput(&input, f);
/*
		if (getInputFloat(&input, "order_parameter", &tmp_value, 0) == KEY_FOUND
			)
			;
	*/
		if(getInputFloat(&input, "RNA_POS_BACK", &tmp_value, 0) == KEY_FOUND) { RNA_POS_BACK = tmp_value; fprintf(stderr,"RNA_POS_BACK  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_POS_STACK", &tmp_value, 0) == KEY_FOUND) { RNA_POS_STACK = tmp_value; fprintf(stderr,"RNA_POS_STACK  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_POS_BASE", &tmp_value, 0) == KEY_FOUND) { RNA_POS_BASE = tmp_value; fprintf(stderr,"RNA_POS_BASE  = %f \n",tmp_value);}

		//NEW STUFF
		if(getInputFloat(&input, "RNA_POS_STACK_3_a1", &tmp_value, 0) == KEY_FOUND) { RNA_POS_STACK_3_a1 = tmp_value; fprintf(stderr,"RNA_POS_STACK_3_a1  = %f \n",tmp_value);}
		if(getInputFloat(&input, "RNA_POS_STACK_3_a2", &tmp_value, 0) == KEY_FOUND) { RNA_POS_STACK_3_a2 = tmp_value; fprintf(stderr,"RNA_POS_STACK_3_a2  = %f \n",tmp_value);}
		if(getInputFloat(&input, "RNA_POS_STACK_5_a1", &tmp_value, 0) == KEY_FOUND) { RNA_POS_STACK_5_a1 = tmp_value; fprintf(stderr,"RNA_POS_STACK_5_a1  = %f \n",tmp_value);}
		if(getInputFloat(&input, "RNA_POS_STACK_5_a2", &tmp_value, 0) == KEY_FOUND) { RNA_POS_STACK_5_a2 = tmp_value; fprintf(stderr,"RNA_POS_STACK_5_a2  = %f \n",tmp_value);}

		if(getInputFloat(&input, "INCLINATION", &tmp_value, 0) == KEY_FOUND) { INCLINATION = tmp_value; fprintf(stderr,"INCLINATION  = %f \n",tmp_value);}
		if(getInputFloat(&input, "BASE_SIZE_R", &tmp_value, 0) == KEY_FOUND) { BASE_SIZE_R = tmp_value; fprintf(stderr,"BASE_SIZE_R  = %f \n",tmp_value);}
		if(getInputFloat(&input, "ROT_SUGAR", &tmp_value, 0) == KEY_FOUND) { ROT_SUGAR = tmp_value; fprintf(stderr,"ROT_SUGAR  = %f \n",tmp_value);}
		if(getInputFloat(&input, "p3_x", &tmp_value, 0) == KEY_FOUND) { p3_x = tmp_value; fprintf(stderr,"p3_x  = %f \n",tmp_value);}
		if(getInputFloat(&input, "p3_y", &tmp_value, 0) == KEY_FOUND) { p3_y = tmp_value; fprintf(stderr,"p3_y  = %f \n",tmp_value);}
		if(getInputFloat(&input, "p3_z", &tmp_value, 0) == KEY_FOUND) { p3_z = tmp_value; fprintf(stderr,"p3_z  = %f \n",tmp_value);}
		if(getInputFloat(&input, "p5_x", &tmp_value, 0) == KEY_FOUND) { p5_x = tmp_value; fprintf(stderr,"p5_x  = %f \n",tmp_value);}
		if(getInputFloat(&input, "p5_y", &tmp_value, 0) == KEY_FOUND) { p5_y = tmp_value; fprintf(stderr,"p5_y  = %f \n",tmp_value);}
		if(getInputFloat(&input, "p5_z", &tmp_value, 0) == KEY_FOUND) { p5_z = tmp_value; fprintf(stderr,"p5_z  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_POS_BACK_a1", &tmp_value, 0) == KEY_FOUND) { RNA_POS_BACK_a1 = tmp_value; fprintf(stderr,"RNA_POS_BACK_a1  = %f \n",tmp_value);}
		if(getInputFloat(&input, "RNA_POS_BACK_a2", &tmp_value, 0) == KEY_FOUND) { RNA_POS_BACK_a2 = tmp_value; fprintf(stderr,"RNA_POS_BACK_a2  = %f \n",tmp_value);}
		if(getInputFloat(&input, "RNA_POS_BACK_a3", &tmp_value, 0) == KEY_FOUND) { RNA_POS_BACK_a3 = tmp_value; fprintf(stderr,"RNA_POS_BACK_a3  = %f \n",tmp_value);}

		// RNA_POS_STACK - RNA_POS_BACK
		if(getInputFloat(&input, "RNA_GAMMA", &tmp_value, 0) == KEY_FOUND) { RNA_GAMMA = tmp_value; fprintf(stderr,"RNA_GAMMA  = %f \n",tmp_value);}


		/*
		 * FENE
		 */
		if(getInputFloat(&input, "RNA_FENE_EPS", &tmp_value, 0) == KEY_FOUND) { RNA_FENE_EPS = tmp_value; fprintf(stderr,"RNA_FENE_EPS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_FENE_R0", &tmp_value, 0) == KEY_FOUND) { RNA_FENE_R0 = tmp_value; fprintf(stderr,"RNA_FENE_R0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_FENE_DELTA", &tmp_value, 0) == KEY_FOUND) { RNA_FENE_DELTA = tmp_value; fprintf(stderr,"RNA_FENE_DELTA  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_FENE_DELTA2", &tmp_value, 0) == KEY_FOUND) { RNA_FENE_DELTA2 = tmp_value; fprintf(stderr,"RNA_FENE_DELTA2  = %f \n",tmp_value);}


		/*
		 * EXCLUDED VOLUME
		 */
		if(getInputFloat(&input, "RNA_EXCL_EPS", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_EPS = tmp_value; fprintf(stderr,"RNA_EXCL_EPS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_S1", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_S1 = tmp_value; fprintf(stderr,"RNA_EXCL_S1  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_S2", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_S2 = tmp_value; fprintf(stderr,"RNA_EXCL_S2  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_S3", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_S3 = tmp_value; fprintf(stderr,"RNA_EXCL_S3  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_S4", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_S4 = tmp_value; fprintf(stderr,"RNA_EXCL_S4  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_R1", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_R1 = tmp_value; fprintf(stderr,"RNA_EXCL_R1  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_R2", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_R2 = tmp_value; fprintf(stderr,"RNA_EXCL_R2  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_R3", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_R3 = tmp_value; fprintf(stderr,"RNA_EXCL_R3  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_R4", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_R4 = tmp_value; fprintf(stderr,"RNA_EXCL_R4  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_B1", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_B1 = tmp_value; fprintf(stderr,"RNA_EXCL_B1  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_B2", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_B2 = tmp_value; fprintf(stderr,"RNA_EXCL_B2  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_B3", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_B3 = tmp_value; fprintf(stderr,"RNA_EXCL_B3  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_B4", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_B4 = tmp_value; fprintf(stderr,"RNA_EXCL_B4  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_RC1", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_RC1 = tmp_value; fprintf(stderr,"RNA_EXCL_RC1  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_RC2", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_RC2 = tmp_value; fprintf(stderr,"RNA_EXCL_RC2  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_RC3", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_RC3 = tmp_value; fprintf(stderr,"RNA_EXCL_RC3  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_EXCL_RC4", &tmp_value, 0) == KEY_FOUND) { RNA_EXCL_RC4 = tmp_value; fprintf(stderr,"RNA_EXCL_RC4  = %f \n",tmp_value);}


		/*
		 * HYDROGEN BONDING
		 */
		// radial part
		if(getInputFloat(&input, "RNA_HYDR_EPS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_EPS = tmp_value; fprintf(stderr,"RNA_HYDR_EPS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_A = tmp_value; fprintf(stderr,"RNA_HYDR_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_RC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_RC = tmp_value; fprintf(stderr,"RNA_HYDR_RC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_R0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_R0 = tmp_value; fprintf(stderr,"RNA_HYDR_R0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_BLOW", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_BLOW = tmp_value; fprintf(stderr,"RNA_HYDR_BLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_BHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_BHIGH = tmp_value; fprintf(stderr,"RNA_HYDR_BHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_RLOW", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_RLOW = tmp_value; fprintf(stderr,"RNA_HYDR_RLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_RHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_RHIGH = tmp_value; fprintf(stderr,"RNA_HYDR_RHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_RCLOW", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_RCLOW = tmp_value; fprintf(stderr,"RNA_HYDR_RCLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_RCHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_RCHIGH = tmp_value; fprintf(stderr,"RNA_HYDR_RCHIGH  = %f \n",tmp_value);}


		// angular part

		if(getInputFloat(&input, "RNA_HYDR_THETA1_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA1_A = tmp_value; fprintf(stderr,"RNA_HYDR_THETA1_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA1_B", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA1_B = tmp_value; fprintf(stderr,"RNA_HYDR_THETA1_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA1_T0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA1_T0 = tmp_value; fprintf(stderr,"RNA_HYDR_THETA1_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA1_TS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA1_TS = tmp_value; fprintf(stderr,"RNA_HYDR_THETA1_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA1_TC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA1_TC = tmp_value; fprintf(stderr,"RNA_HYDR_THETA1_TC  = %f \n",tmp_value);}


		if(getInputFloat(&input, "RNA_HYDR_THETA2_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA2_A = tmp_value; fprintf(stderr,"RNA_HYDR_THETA2_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA2_B", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA2_B = tmp_value; fprintf(stderr,"RNA_HYDR_THETA2_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA2_T0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA2_T0 = tmp_value; fprintf(stderr,"RNA_HYDR_THETA2_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA2_TS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA2_TS = tmp_value; fprintf(stderr,"RNA_HYDR_THETA2_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA2_TC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA2_TC = tmp_value; fprintf(stderr,"RNA_HYDR_THETA2_TC  = %f \n",tmp_value);}


		if(getInputFloat(&input, "RNA_HYDR_THETA3_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA3_A = tmp_value; fprintf(stderr,"RNA_HYDR_THETA3_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA3_B", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA3_B = tmp_value; fprintf(stderr,"RNA_HYDR_THETA3_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA3_T0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA3_T0 = tmp_value; fprintf(stderr,"RNA_HYDR_THETA3_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA3_TS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA3_TS = tmp_value; fprintf(stderr,"RNA_HYDR_THETA3_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA3_TC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA3_TC = tmp_value; fprintf(stderr,"RNA_HYDR_THETA3_TC  = %f \n",tmp_value);}


		if(getInputFloat(&input, "RNA_HYDR_THETA4_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA4_A = tmp_value; fprintf(stderr,"RNA_HYDR_THETA4_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA4_B", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA4_B = tmp_value; fprintf(stderr,"RNA_HYDR_THETA4_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA4_T0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA4_T0 = tmp_value; fprintf(stderr,"RNA_HYDR_THETA4_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA4_TS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA4_TS = tmp_value; fprintf(stderr,"RNA_HYDR_THETA4_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA4_TC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA4_TC = tmp_value; fprintf(stderr,"RNA_HYDR_THETA4_TC  = %f \n",tmp_value);}


		if(getInputFloat(&input, "RNA_HYDR_THETA7_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA7_A = tmp_value; fprintf(stderr,"RNA_HYDR_THETA7_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA7_B", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA7_B = tmp_value; fprintf(stderr,"RNA_HYDR_THETA7_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA7_T0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA7_T0 = tmp_value; fprintf(stderr,"RNA_HYDR_THETA7_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA7_TS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA7_TS = tmp_value; fprintf(stderr,"RNA_HYDR_THETA7_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA7_TC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA7_TC = tmp_value; fprintf(stderr,"RNA_HYDR_THETA7_TC  = %f \n",tmp_value);}


		if(getInputFloat(&input, "RNA_HYDR_THETA8_A", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA8_A = tmp_value; fprintf(stderr,"RNA_HYDR_THETA8_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA8_B", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA8_B = tmp_value; fprintf(stderr,"RNA_HYDR_THETA8_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA8_T0", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA8_T0 = tmp_value; fprintf(stderr,"RNA_HYDR_THETA8_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA8_TS", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA8_TS = tmp_value; fprintf(stderr,"RNA_HYDR_THETA8_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_HYDR_THETA8_TC", &tmp_value, 0) == KEY_FOUND) { RNA_HYDR_THETA8_TC = tmp_value; fprintf(stderr,"RNA_HYDR_THETA8_TC  = %f \n",tmp_value);}


		/*
		 * STACKING
		 */
		// radial part
		if(getInputFloat(&input, "RNA_STCK_BASE_EPS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_BASE_EPS = tmp_value; fprintf(stderr,"RNA_STCK_BASE_EPS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_FACT_EPS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_FACT_EPS = tmp_value; fprintf(stderr,"RNA_STCK_FACT_EPS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_A", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_A = tmp_value; fprintf(stderr,"RNA_STCK_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_RC", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_RC = tmp_value; fprintf(stderr,"RNA_STCK_RC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_R0", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_R0 = tmp_value; fprintf(stderr,"RNA_STCK_R0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_BLOW", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_BLOW = tmp_value; fprintf(stderr,"RNA_STCK_BLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_BHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_BHIGH = tmp_value; fprintf(stderr,"RNA_STCK_BHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_RLOW", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_RLOW = tmp_value; fprintf(stderr,"RNA_STCK_RLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_RHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_RHIGH = tmp_value; fprintf(stderr,"RNA_STCK_RHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_RCLOW", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_RCLOW = tmp_value; fprintf(stderr,"RNA_STCK_RCLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_RCHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_RCHIGH = tmp_value; fprintf(stderr,"RNA_STCK_RCHIGH  = %f \n",tmp_value);}


		// angular part, theta
		if(getInputFloat(&input, "RNA_STCK_THETA4_A", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA4_A = tmp_value; fprintf(stderr,"RNA_STCK_THETA4_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA4_B", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA4_B = tmp_value; fprintf(stderr,"RNA_STCK_THETA4_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA4_T0", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA4_T0 = tmp_value; fprintf(stderr,"RNA_STCK_THETA4_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA4_TS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA4_TS = tmp_value; fprintf(stderr,"RNA_STCK_THETA4_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA4_TC", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA4_TC = tmp_value; fprintf(stderr,"RNA_STCK_THETA4_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA5_A", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA5_A = tmp_value; fprintf(stderr,"RNA_STCK_THETA5_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA5_B", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA5_B = tmp_value; fprintf(stderr,"RNA_STCK_THETA5_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA5_T0", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA5_T0 = tmp_value; fprintf(stderr,"RNA_STCK_THETA5_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA5_TS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA5_TS = tmp_value; fprintf(stderr,"RNA_STCK_THETA5_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA5_TC", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA5_TC = tmp_value; fprintf(stderr,"RNA_STCK_THETA5_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA6_A", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA6_A = tmp_value; fprintf(stderr,"RNA_STCK_THETA6_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA6_B", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA6_B = tmp_value; fprintf(stderr,"RNA_STCK_THETA6_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA6_T0", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA6_T0 = tmp_value; fprintf(stderr,"RNA_STCK_THETA6_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA6_TS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA6_TS = tmp_value; fprintf(stderr,"RNA_STCK_THETA6_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_THETA6_TC", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_THETA6_TC = tmp_value; fprintf(stderr,"RNA_STCK_THETA6_TC  = %f \n",tmp_value);}


		//ANGULAR MODULATION B1, B2
		if(getInputFloat(&input, "STCK_THETAB1_A", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB1_A = tmp_value; fprintf(stderr,"STCK_THETAB1_A  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB1_B", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB1_B = tmp_value; fprintf(stderr,"STCK_THETAB1_B  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB1_T0", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB1_T0 = tmp_value; fprintf(stderr,"STCK_THETAB1_T0  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB1_TS", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB1_TS = tmp_value; fprintf(stderr,"STCK_THETAB1_TS  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB1_TC", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB1_TC = tmp_value; fprintf(stderr,"STCK_THETAB1_TC  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB2_A", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB2_A = tmp_value; fprintf(stderr,"STCK_THETAB2_A  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB2_B", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB2_B = tmp_value; fprintf(stderr,"STCK_THETAB2_B  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB2_T0", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB2_T0 = tmp_value; fprintf(stderr,"STCK_THETAB2_T0  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB2_TS", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB2_TS = tmp_value; fprintf(stderr,"STCK_THETAB2_TS  = %f \n",tmp_value);}

				if(getInputFloat(&input, "STCK_THETAB2_TC", &tmp_value, 0) == KEY_FOUND) { STCK_THETAB2_TC = tmp_value; fprintf(stderr,"STCK_THETAB2_TC  = %f \n",tmp_value);}





		// angular part, phi
		if(getInputFloat(&input, "RNA_STCK_PHI1_A", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI1_A = tmp_value; fprintf(stderr,"RNA_STCK_PHI1_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI1_B", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI1_B = tmp_value; fprintf(stderr,"RNA_STCK_PHI1_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI1_XC", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI1_XC = tmp_value; fprintf(stderr,"RNA_STCK_PHI1_XC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI1_XS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI1_XS = tmp_value; fprintf(stderr,"RNA_STCK_PHI1_XS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI2_A", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI2_A = tmp_value; fprintf(stderr,"RNA_STCK_PHI2_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI2_B", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI2_B = tmp_value; fprintf(stderr,"RNA_STCK_PHI2_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI2_XC", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI2_XC = tmp_value; fprintf(stderr,"RNA_STCK_PHI2_XC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_STCK_PHI2_XS", &tmp_value, 0) == KEY_FOUND) { RNA_STCK_PHI2_XS = tmp_value; fprintf(stderr,"RNA_STCK_PHI2_XS  = %f \n",tmp_value);}


		/*
		 *  CROSS STACKING
		 */
		// radial part
		if(getInputFloat(&input, "RNA_CRST_R0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_R0 = tmp_value; fprintf(stderr,"RNA_CRST_R0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_RC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_RC = tmp_value; fprintf(stderr,"RNA_CRST_RC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_K", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_K = tmp_value; fprintf(stderr,"RNA_CRST_K  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_BLOW", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_BLOW = tmp_value; fprintf(stderr,"RNA_CRST_BLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_RLOW", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_RLOW = tmp_value; fprintf(stderr,"RNA_CRST_RLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_RCLOW", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_RCLOW = tmp_value; fprintf(stderr,"RNA_CRST_RCLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_BHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_BHIGH = tmp_value; fprintf(stderr,"RNA_CRST_BHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_RHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_RHIGH = tmp_value; fprintf(stderr,"RNA_CRST_RHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_RCHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_RCHIGH = tmp_value; fprintf(stderr,"RNA_CRST_RCHIGH  = %f \n",tmp_value);}

		// angular part;
		if(getInputFloat(&input, "RNA_CRST_THETA1_A", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA1_A = tmp_value; fprintf(stderr,"RNA_CRST_THETA1_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA1_B", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA1_B = tmp_value; fprintf(stderr,"RNA_CRST_THETA1_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA1_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA1_T0 = tmp_value; fprintf(stderr,"RNA_CRST_THETA1_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA1_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA1_TS = tmp_value; fprintf(stderr,"RNA_CRST_THETA1_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA1_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA1_TC = tmp_value; fprintf(stderr,"RNA_CRST_THETA1_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA2_A", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA2_A = tmp_value; fprintf(stderr,"RNA_CRST_THETA2_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA2_B", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA2_B = tmp_value; fprintf(stderr,"RNA_CRST_THETA2_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA2_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA2_T0 = tmp_value; fprintf(stderr,"RNA_CRST_THETA2_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA2_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA2_TS = tmp_value; fprintf(stderr,"RNA_CRST_THETA2_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA2_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA2_TC = tmp_value; fprintf(stderr,"RNA_CRST_THETA2_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA3_A", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA3_A = tmp_value; fprintf(stderr,"RNA_CRST_THETA3_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA3_B", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA3_B = tmp_value; fprintf(stderr,"RNA_CRST_THETA3_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA3_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA3_T0 = tmp_value; fprintf(stderr,"RNA_CRST_THETA3_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA3_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA3_TS = tmp_value; fprintf(stderr,"RNA_CRST_THETA3_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA3_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA3_TC = tmp_value; fprintf(stderr,"RNA_CRST_THETA3_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA4_A", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA4_A = tmp_value; fprintf(stderr,"RNA_CRST_THETA4_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA4_B", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA4_B = tmp_value; fprintf(stderr,"RNA_CRST_THETA4_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA4_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA4_T0 = tmp_value; fprintf(stderr,"RNA_CRST_THETA4_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA4_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA4_TS = tmp_value; fprintf(stderr,"RNA_CRST_THETA4_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA4_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA4_TC = tmp_value; fprintf(stderr,"RNA_CRST_THETA4_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA7_A", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA7_A = tmp_value; fprintf(stderr,"RNA_CRST_THETA7_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA7_B", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA7_B = tmp_value; fprintf(stderr,"RNA_CRST_THETA7_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA7_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA7_T0 = tmp_value; fprintf(stderr,"RNA_CRST_THETA7_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA7_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA7_TS = tmp_value; fprintf(stderr,"RNA_CRST_THETA7_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA7_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA7_TC = tmp_value; fprintf(stderr,"RNA_CRST_THETA7_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA8_A", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA8_A = tmp_value; fprintf(stderr,"RNA_CRST_THETA8_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA8_B", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA8_B = tmp_value; fprintf(stderr,"RNA_CRST_THETA8_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA8_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA8_T0 = tmp_value; fprintf(stderr,"RNA_CRST_THETA8_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA8_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA8_TS = tmp_value; fprintf(stderr,"RNA_CRST_THETA8_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CRST_THETA8_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CRST_THETA8_TC = tmp_value; fprintf(stderr,"RNA_CRST_THETA8_TC  = %f \n",tmp_value);}


		/*
		 *  COAXIAL STACKING
		 */
		// radial part
		if(getInputFloat(&input, "RNA_CXST_R0", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_R0 = tmp_value; fprintf(stderr,"RNA_CXST_R0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_RC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_RC = tmp_value; fprintf(stderr,"RNA_CXST_RC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_K", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_K = tmp_value; fprintf(stderr,"RNA_CXST_K  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_BLOW", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_BLOW = tmp_value; fprintf(stderr,"RNA_CXST_BLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_RLOW", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_RLOW = tmp_value; fprintf(stderr,"RNA_CXST_RLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_RCLOW", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_RCLOW = tmp_value; fprintf(stderr,"RNA_CXST_RCLOW  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_BHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_BHIGH = tmp_value; fprintf(stderr,"RNA_CXST_BHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_RHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_RHIGH = tmp_value; fprintf(stderr,"RNA_CXST_RHIGH  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_RCHIGH", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_RCHIGH = tmp_value; fprintf(stderr,"RNA_CXST_RCHIGH  = %f \n",tmp_value);}

		// angular part;

		if(getInputFloat(&input, "RNA_CXST_THETA1_A", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA1_A = tmp_value; fprintf(stderr,"RNA_CXST_THETA1_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA1_B", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA1_B = tmp_value; fprintf(stderr,"RNA_CXST_THETA1_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA1_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA1_T0 = tmp_value; fprintf(stderr,"RNA_CXST_THETA1_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA1_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA1_TS = tmp_value; fprintf(stderr,"RNA_CXST_THETA1_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA1_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA1_TC = tmp_value; fprintf(stderr,"RNA_CXST_THETA1_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA4_A", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA4_A = tmp_value; fprintf(stderr,"RNA_CXST_THETA4_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA4_B", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA4_B = tmp_value; fprintf(stderr,"RNA_CXST_THETA4_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA4_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA4_T0 = tmp_value; fprintf(stderr,"RNA_CXST_THETA4_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA4_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA4_TS = tmp_value; fprintf(stderr,"RNA_CXST_THETA4_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA4_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA4_TC = tmp_value; fprintf(stderr,"RNA_CXST_THETA4_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA5_A", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA5_A = tmp_value; fprintf(stderr,"RNA_CXST_THETA5_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA5_B", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA5_B = tmp_value; fprintf(stderr,"RNA_CXST_THETA5_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA5_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA5_T0 = tmp_value; fprintf(stderr,"RNA_CXST_THETA5_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA5_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA5_TS = tmp_value; fprintf(stderr,"RNA_CXST_THETA5_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA5_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA5_TC = tmp_value; fprintf(stderr,"RNA_CXST_THETA5_TC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA6_A", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA6_A = tmp_value; fprintf(stderr,"RNA_CXST_THETA6_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA6_B", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA6_B = tmp_value; fprintf(stderr,"RNA_CXST_THETA6_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA6_T0", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA6_T0 = tmp_value; fprintf(stderr,"RNA_CXST_THETA6_T0  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA6_TS", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA6_TS = tmp_value; fprintf(stderr,"RNA_CXST_THETA6_TS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_THETA6_TC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_THETA6_TC = tmp_value; fprintf(stderr,"RNA_CXST_THETA6_TC  = %f \n",tmp_value);}

		// angular part, phi3 and phi4

		if(getInputFloat(&input, "RNA_CXST_PHI3_A", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI3_A = tmp_value; fprintf(stderr,"RNA_CXST_PHI3_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI3_B", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI3_B = tmp_value; fprintf(stderr,"RNA_CXST_PHI3_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI3_XC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI3_XC = tmp_value; fprintf(stderr,"RNA_CXST_PHI3_XC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI3_XS", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI3_XS = tmp_value; fprintf(stderr,"RNA_CXST_PHI3_XS  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI4_A", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI4_A = tmp_value; fprintf(stderr,"RNA_CXST_PHI4_A  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI4_B", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI4_B = tmp_value; fprintf(stderr,"RNA_CXST_PHI4_B  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI4_XC", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI4_XC = tmp_value; fprintf(stderr,"RNA_CXST_PHI4_XC  = %f \n",tmp_value);}

		if(getInputFloat(&input, "RNA_CXST_PHI4_XS", &tmp_value, 0) == KEY_FOUND) { RNA_CXST_PHI4_XS = tmp_value; fprintf(stderr,"RNA_CXST_PHI4_XS  = %f \n",tmp_value);}


		fclose(f);
	}

};
/*
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
*/
#endif /* MODEL_H_ */
